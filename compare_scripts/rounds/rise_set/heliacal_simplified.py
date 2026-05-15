#!/usr/bin/env python3
"""Round 66: Heliacal Events (simplified — heliacal calcs are very slow)"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0

DATM = (1013.25, 15.0, 40.0, 0.0)
DOBS = (36.0, 1.0, 1, 0, 0, 0)
geopos = (12.50, 41.90, 0.0)  # Rome
jd_start = 2451545.0

print("=" * 70)
print("ROUND 66: Heliacal Events (Simplified)")
print("=" * 70)

# P1: Rising for visible planets
print("\n=== P1: Heliacal Rising ===")
for planet_name in ["Mercury", "Venus", "Mars", "Jupiter", "Saturn"]:
    try:
        se_r = swe.heliacal_ut(jd_start, geopos, DATM, DOBS, planet_name, 1, 2)
        le_r = ephem.heliacal_ut(jd_start, geopos, DATM, DOBS, planet_name, 1, 2)
        diff = abs(se_r[0] - le_r[0])
        if diff < 2.0:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL {planet_name}: SE={se_r[0]:.4f} LE={le_r[0]:.4f} diff={diff:.2f}d"
            )
    except Exception as e:
        errors += 1
        print(f"  ERR {planet_name}: {str(e)[:80]}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Setting for visible planets
print("\n=== P2: Heliacal Setting ===")
for planet_name in ["Mercury", "Venus", "Mars", "Jupiter", "Saturn"]:
    try:
        se_r = swe.heliacal_ut(jd_start, geopos, DATM, DOBS, planet_name, 2, 2)
        le_r = ephem.heliacal_ut(jd_start, geopos, DATM, DOBS, planet_name, 2, 2)
        diff = abs(se_r[0] - le_r[0])
        if diff < 2.0:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL {planet_name}: SE={se_r[0]:.4f} LE={le_r[0]:.4f} diff={diff:.2f}d"
            )
    except Exception as e:
        errors += 1
        print(f"  ERR {planet_name}: {str(e)[:80]}")
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Evening First / Morning Last for inner planets
print("\n=== P3: Evening First / Morning Last ===")
for planet_name in ["Mercury", "Venus"]:
    for evt in [3, 4]:
        try:
            se_r = swe.heliacal_ut(jd_start, geopos, DATM, DOBS, planet_name, evt, 2)
            le_r = ephem.heliacal_ut(
                jd_start, geopos, DATM, DOBS, planet_name, evt, 2
            )
            diff = abs(se_r[0] - le_r[0])
            if diff < 2.0:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL {planet_name} evt={evt}: SE={se_r[0]:.4f} LE={le_r[0]:.4f} diff={diff:.2f}d"
                )
        except Exception as e:
            errors += 1
            print(f"  ERR {planet_name} evt={evt}: {str(e)[:80]}")
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Stars
print("\n=== P4: Fixed star rising ===")
for star in ["Sirius", "Aldebaran", "Regulus", "Spica", "Antares"]:
    try:
        se_r = swe.heliacal_ut(jd_start, geopos, DATM, DOBS, star, 1, 2)
        le_r = ephem.heliacal_ut(jd_start, geopos, DATM, DOBS, star, 1, 2)
        diff = abs(se_r[0] - le_r[0])
        if diff < 3.0:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL {star}: SE={se_r[0]:.4f} LE={le_r[0]:.4f} diff={diff:.2f}d")
    except Exception as e:
        errors += 1
        print(f"  ERR {star}: {str(e)[:80]}")
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Return value consistency
print("\n=== P5: Return value consistency ===")
for planet_name in ["Venus", "Jupiter"]:
    try:
        le_r = ephem.heliacal_ut(jd_start, geopos, DATM, DOBS, planet_name, 1, 2)
        if le_r[0] > 0 and le_r[1] > 0 and le_r[2] > 0:
            if le_r[0] <= le_r[1] + 0.01 and le_r[1] <= le_r[2] + 0.01:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL {planet_name}: start={le_r[0]:.4f} opt={le_r[1]:.4f} end={le_r[2]:.4f}"
                )
        else:
            passed += 1
    except:
        errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 66 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
