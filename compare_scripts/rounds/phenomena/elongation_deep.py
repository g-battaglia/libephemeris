#!/usr/bin/env python3
"""Round 67: Planetary Elongation & Phenomena Deep Sweep"""

from __future__ import annotations
import sys, os, math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0

BODIES = list(range(1, 10))  # Moon-Pluto (skip Sun)
NAMES = [
    "_",
    "Moon",
    "Mercury",
    "Venus",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "Pluto",
]
FLAGS = 256
DATES = [2451545.0 + i * 450 for i in range(20)]

print("=" * 70)
print("ROUND 67: Planetary Elongation & Phenomena Deep Sweep")
print("=" * 70)

# P1: Elongation
print("\n=== P1: Elongation ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.pheno_ut(jd, body, FLAGS)  # flat tuple
            le = ephem.pheno_ut(jd, body, FLAGS)  # (tuple, flag)
            diff = abs(se[2] - le[0][2])
            if diff < 0.01:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P1 {NAMES[body]} jd={jd:.0f}: SE={se[2]:.6f} LE={le[0][2]:.6f} d={diff:.6f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Phase angle
print("\n=== P2: Phase angle ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.pheno_ut(jd, body, FLAGS)
            le = ephem.pheno_ut(jd, body, FLAGS)
            diff = abs(se[0] - le[0][0])
            if diff < 0.5:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P2 {NAMES[body]} jd={jd:.0f}: SE={se[0]:.4f} LE={le[0][0]:.4f} d={diff:.4f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Phase (illuminated fraction)
print("\n=== P3: Phase (illuminated fraction) ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.pheno_ut(jd, body, FLAGS)
            le = ephem.pheno_ut(jd, body, FLAGS)
            diff = abs(se[1] - le[0][1])
            if diff < 0.005:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P3 {NAMES[body]} jd={jd:.0f}: SE={se[1]:.6f} LE={le[0][1]:.6f} d={diff:.6f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Apparent diameter
print("\n=== P4: Apparent diameter ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.pheno_ut(jd, body, FLAGS)
            le = ephem.pheno_ut(jd, body, FLAGS)
            if se[3] > 0 and le[0][3] > 0:
                ratio = le[0][3] / se[3]
                if 0.95 < ratio < 1.05:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P4 {NAMES[body]} jd={jd:.0f}: SE={se[3]:.6f} LE={le[0][3]:.6f} r={ratio:.4f}"
                    )
            else:
                passed += 1
        except Exception as e:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Magnitude
print("\n=== P5: Magnitude ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.pheno_ut(jd, body, FLAGS)
            le = ephem.pheno_ut(jd, body, FLAGS)
            diff = abs(se[4] - le[0][4])
            if diff < 0.5:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P5 {NAMES[body]} jd={jd:.0f}: SE={se[4]:.4f} LE={le[0][4]:.4f} d={diff:.4f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Manual elongation cross-check
print("\n=== P6: Manual elongation cross-check ===")
for jd in DATES[:10]:
    for body in [1, 2, 3, 4, 5]:
        try:
            sun = ephem.calc_ut(jd, 0, FLAGS)[0]
            bdy = ephem.calc_ut(jd, body, FLAGS)[0]
            dlon = math.radians(bdy[0] - sun[0])
            cos_e = math.sin(math.radians(sun[1])) * math.sin(
                math.radians(bdy[1])
            ) + math.cos(math.radians(sun[1])) * math.cos(
                math.radians(bdy[1])
            ) * math.cos(dlon)
            manual = math.degrees(math.acos(max(-1.0, min(1.0, cos_e))))
            le = ephem.pheno_ut(jd, body, FLAGS)
            diff = abs(manual - le[0][2])
            if diff < 0.1:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P6 {NAMES[body]} jd={jd:.0f}: manual={manual:.4f} pheno={le[0][2]:.4f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 67 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
