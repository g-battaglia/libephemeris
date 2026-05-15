#!/usr/bin/env python3
"""Round 76: Declination Parallels & Out-of-Bounds Planets"""

from __future__ import annotations
import sys, os, math
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256 | 2048  # FLG_SPEED | FLG_EQUATORIAL

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
DATES = [2451545.0 + i * 30 for i in range(240)]  # 20 years monthly

print("=" * 70)
print("ROUND 76: Declination Parallels & Out-of-Bounds Planets")
print("=" * 70)

# P1: Declination comparison for all planets
print("\n=== P1: Declination all planets (20-year monthly) ===")
for body in range(10):
    body_pass = body_fail = 0
    for jd in DATES:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            se_dec = se[0][1]
            le_dec = le[0][1]
            diff = abs(se_dec - le_dec) * 3600
            if diff < 1.0:
                passed += 1
                body_pass += 1
            else:
                failed += 1
                body_fail += 1
        except:
            errors += 1
    if body_fail > 0:
        print(f"  {NAMES[body]}: {body_pass}/{body_pass + body_fail}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Out-of-bounds detection (dec > 23.44° or < -23.44°)
print("\n=== P2: Out-of-bounds detection ===")
obliquity = 23.44
oob_count = 0
oob_agree = 0
for body in range(10):
    for jd in DATES:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            se_oob = abs(se[0][1]) > obliquity
            le_oob = abs(le[0][1]) > obliquity
            if se_oob or le_oob:
                oob_count += 1
                if se_oob == le_oob:
                    passed += 1
                    oob_agree += 1
                else:
                    failed += 1
            else:
                passed += 1
        except:
            errors += 1
print(f"  OOB events found: {oob_count}, agreement: {oob_agree}")
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Right Ascension comparison
print("\n=== P3: Right Ascension all planets ===")
for body in range(10):
    for jd in DATES[:120]:  # 10 years
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
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

# P4: Declination speed
print("\n=== P4: Declination speed ===")
for body in range(10):
    for jd in DATES[:60]:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][4] - le[0][4])
            if diff < 0.001:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Moon maximum declination (should reach ~28.6° at max)
print("\n=== P5: Moon max declination check ===")
max_se_dec = 0
max_le_dec = 0
for i in range(6935):  # ~19 years (one nodal cycle) daily
    jd = 2451545.0 + i
    try:
        se = swe.calc_ut(jd, 1, FLAGS)
        le = ephem.calc_ut(jd, 1, FLAGS)
        max_se_dec = max(max_se_dec, abs(se[0][1]))
        max_le_dec = max(max_le_dec, abs(le[0][1]))
    except:
        pass
diff = abs(max_se_dec - max_le_dec)
if diff < 0.01:
    passed += 1
    print(f"  Moon max dec: SE={max_se_dec:.4f}° LE={max_le_dec:.4f}° diff={diff:.6f}°")
else:
    failed += 1
    print(
        f"  FAIL Moon max dec: SE={max_se_dec:.4f}° LE={max_le_dec:.4f}° diff={diff:.4f}°"
    )
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Venus OOB verification (Venus can go to ~27°)
print("\n=== P6: Venus max declination ===")
max_venus = 0
for i in range(2920):  # 8 years
    jd = 2451545.0 + i
    try:
        le = ephem.calc_ut(jd, 3, FLAGS)
        max_venus = max(max_venus, abs(le[0][1]))
    except:
        pass
# Venus should reach at least 24° at some point
if max_venus > 24.0:
    passed += 1
    print(f"  Venus max dec: {max_venus:.4f}° (OOB confirmed)")
else:
    failed += 1
    print(f"  FAIL Venus max dec: {max_venus:.4f}° (expected > 24°)")
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 76 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
