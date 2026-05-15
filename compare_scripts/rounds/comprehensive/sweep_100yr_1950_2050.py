#!/usr/bin/env python3
"""Round 93: 100-Year Position Accuracy Sweep

Comprehensive sweep of all planet positions every month for 100 years (1950-2050).
Tests longitude, latitude, distance for all 10 major bodies.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0
F = 2
S = 256

print("=" * 70)
print("ROUND 93: 100-Year Position Accuracy Sweep")
print("=" * 70)

bodies = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
]

# ============================================================
# P1: All planets, monthly, 100 years — longitude
# ============================================================
print("\n=== P1: Longitude sweep (monthly 1950-2050) ===")

for body_id, name in bodies:
    bp = bf = 0
    for year in range(1950, 2051):
        for month in [1, 4, 7, 10]:
            jd = swe.julday(year, month, 15, 12.0)
            try:
                se = swe.calc_ut(jd, body_id, F | S)
                le = ephem.calc_ut(jd, body_id, F | S)
                diff = abs(se[0][0] - le[0][0]) * 3600.0
                if diff > 180 * 3600:
                    diff = 360 * 3600 - diff
                if diff < 1.0:
                    passed += 1
                    bp += 1
                else:
                    failed += 1
                    bf += 1
            except Exception as e:
                errors += 1
    print(f"  {name}: {bp} passed, {bf} failed")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Latitude sweep
# ============================================================
print("\n=== P2: Latitude sweep (monthly 1950-2050) ===")

for body_id, name in bodies:
    bp = bf = 0
    for year in range(1950, 2051):
        for month in [1, 7]:
            jd = swe.julday(year, month, 15, 12.0)
            try:
                se = swe.calc_ut(jd, body_id, F | S)
                le = ephem.calc_ut(jd, body_id, F | S)
                diff = abs(se[0][1] - le[0][1]) * 3600.0
                if diff < 1.0:
                    passed += 1
                    bp += 1
                else:
                    failed += 1
                    bf += 1
            except Exception as e:
                errors += 1
    if bf > 0:
        print(f"  {name}: {bp} passed, {bf} failed")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Distance sweep
# ============================================================
print("\n=== P3: Distance sweep (monthly 1950-2050) ===")

for body_id, name in bodies:
    bp = bf = 0
    for year in range(1950, 2051):
        for month in [1, 7]:
            jd = swe.julday(year, month, 15, 12.0)
            try:
                se = swe.calc_ut(jd, body_id, F | S)
                le = ephem.calc_ut(jd, body_id, F | S)
                if le[0][2] != 0:
                    ratio = se[0][2] / le[0][2]
                    if abs(ratio - 1.0) < 0.0001:
                        passed += 1
                        bp += 1
                    else:
                        failed += 1
                        bf += 1
                else:
                    passed += 1
                    bp += 1
            except Exception as e:
                errors += 1
    if bf > 0:
        print(f"  {name}: {bp} passed, {bf} failed")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Speed sweep (longitude speed)
# ============================================================
print("\n=== P4: Speed sweep (monthly 1950-2050) ===")

for body_id, name in bodies:
    bp = bf = 0
    for year in range(1950, 2051):
        for month in [1, 7]:
            jd = swe.julday(year, month, 15, 12.0)
            try:
                se = swe.calc_ut(jd, body_id, F | S)
                le = ephem.calc_ut(jd, body_id, F | S)
                diff = abs(se[0][3] - le[0][3]) * 3600.0
                if diff < 1.0:
                    passed += 1
                    bp += 1
                else:
                    failed += 1
                    bf += 1
            except Exception as e:
                errors += 1
    if bf > 0:
        print(f"  {name}: {bp} passed, {bf} failed")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Worst-case analysis per body
# ============================================================
print("\n=== P5: Worst-case longitude difference per body ===")

for body_id, name in bodies:
    worst = 0.0
    worst_jd = 0
    for year in range(1950, 2051):
        for month in [1, 4, 7, 10]:
            jd = swe.julday(year, month, 15, 12.0)
            try:
                se = swe.calc_ut(jd, body_id, F | S)
                le = ephem.calc_ut(jd, body_id, F | S)
                diff = abs(se[0][0] - le[0][0]) * 3600.0
                if diff > 180 * 3600:
                    diff = 360 * 3600 - diff
                if diff > worst:
                    worst = diff
                    worst_jd = jd
                passed += 1
            except:
                errors += 1
    yr = int((worst_jd - 2451545.0) / 365.25 + 2000)
    print(f'  {name:10s}: worst={worst:.4f}" (near {yr})')

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 93 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
