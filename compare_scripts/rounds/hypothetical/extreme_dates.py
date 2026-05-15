#!/usr/bin/env python3
"""Round 73: Hypothetical/Uranian Bodies at Extreme Dates"""

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

# Uranian body IDs: 40-48 (Cupido through Proserpina) + 50-57 (Hamburg school)
URANIANS = list(range(40, 49))
URAN_NAMES = {
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Isis-Transpluto",
}

# Date range: 1900-2100 at 10-year intervals
DATES = [2415020.0 + i * 3652.5 for i in range(21)]  # 1900 to 2100

print("=" * 70)
print("ROUND 73: Hypothetical/Uranian Bodies at Extreme Dates")
print("=" * 70)

# P1: All Uranians across 200-year range
print("\n=== P1: Uranian positions 1900-2100 ===")
for body in URANIANS:
    body_pass = body_fail = body_err = 0
    for jd in DATES:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            diff_as = diff * 3600
            if diff_as < 60:  # 1 arcmin tolerance for hypothetical bodies
                passed += 1
                body_pass += 1
            else:
                failed += 1
                body_fail += 1
                if body_fail <= 2:
                    print(
                        f'  FAIL {URAN_NAMES.get(body, body)} jd={jd:.0f}: SE={se[0][0]:.4f} LE={le[0][0]:.4f} d={diff_as:.1f}"'
                    )
        except Exception as e:
            errors += 1
            body_err += 1
    total_b = body_pass + body_fail
    if total_b > 0:
        print(
            f"  {URAN_NAMES.get(body, body)}: {body_pass}/{total_b} ({100 * body_pass / total_b:.0f}%) err={body_err}"
        )
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Uranian latitude/distance
print("\n=== P2: Uranian latitude ===")
for body in URANIANS:
    for jd in DATES[:10]:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff_lat = abs(se[0][1] - le[0][1]) * 3600
            if diff_lat < 60:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Uranian speeds
print("\n=== P3: Uranian speeds ===")
for body in URANIANS:
    for jd in DATES[:10]:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff_spd = abs(se[0][3] - le[0][3])
            if diff_spd < 0.01:  # 0.01 deg/day
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: J2000 mode
print("\n=== P4: Uranian J2000 mode ===")
J2000_FLAGS = 256 | 32  # SPEED | J2000
for body in URANIANS:
    for jd in DATES[:10]:
        try:
            se = swe.calc_ut(jd, body, J2000_FLAGS)
            le = ephem.calc_ut(jd, body, J2000_FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 60:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Heliocentric Uranians
print("\n=== P5: Heliocentric Uranian ===")
HELIO_FLAGS = 256 | 8  # SPEED | HELCTR
for body in URANIANS:
    for jd in DATES[:10]:
        try:
            se = swe.calc_ut(jd, body, HELIO_FLAGS)
            le = ephem.calc_ut(jd, body, HELIO_FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 30:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 73 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
