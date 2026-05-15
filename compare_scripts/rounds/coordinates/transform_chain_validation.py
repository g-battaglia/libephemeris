#!/usr/bin/env python3
"""Round 72: Coordinate Transformation Chain Validation"""

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

DATES = [2451545.0 + i * 365.25 for i in range(20)]
BODIES = list(range(10))
NAMES = [
    "Sun",
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

print("=" * 70)
print("ROUND 72: Coordinate Transformation Chain Validation")
print("=" * 70)

# P1: cotrans ecliptic→equatorial→ecliptic round-trip
print("\n=== P1: cotrans round-trip ===")
for jd in DATES[:10]:
    obl = ephem.calc_ut(jd, -1, 0)[0][1]
    for body in BODIES:
        try:
            pos = ephem.calc_ut(jd, body, 256)[0]
            lon, lat = pos[0], pos[1]
            eq = ephem.cotrans((lon, lat, 1.0), -obl)
            ecl = ephem.cotrans((eq[0], eq[1], 1.0), obl)
            diff_lon = abs(ecl[0] - lon)
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lat = abs(ecl[1] - lat)
            if diff_lon < 1e-8 and diff_lat < 1e-8:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(
                        f"  FAIL P1 {NAMES[body]} jd={jd:.0f}: dlon={diff_lon:.2e} dlat={diff_lat:.2e}"
                    )
        except Exception as e:
            errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: FLG_EQUATORIAL vs manual cotrans
print("\n=== P2: FLG_EQUATORIAL vs manual cotrans ===")
for jd in DATES:
    obl = ephem.calc_ut(jd, -1, 0)[0][1]
    for body in BODIES:
        try:
            ecl = ephem.calc_ut(jd, body, 256)[0]
            eq_flag = ephem.calc_ut(jd, body, 256 | 2048)[0]
            eq_manual = ephem.cotrans((ecl[0], ecl[1], 1.0), -obl)
            diff_ra = abs(eq_flag[0] - eq_manual[0])
            if diff_ra > 180:
                diff_ra = 360 - diff_ra
            diff_dec = abs(eq_flag[1] - eq_manual[1])
            if diff_ra * 3600 < 0.1 and diff_dec * 3600 < 0.1:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(
                        f'  FAIL P2 {NAMES[body]} jd={jd:.0f}: dra={diff_ra * 3600:.3f}" ddec={diff_dec * 3600:.3f}"'
                    )
        except Exception as e:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: SE vs LE cotrans comparison
print("\n=== P3: SE vs LE cotrans ===")
test_coords = [
    (0, 0),
    (45, 10),
    (90, 23),
    (135, -15),
    (180, 5),
    (225, -20),
    (270, 10),
    (315, 0),
]
for lon, lat in test_coords:
    for obl in [23.0, 23.4393, 23.5, 24.0]:
        try:
            se_eq = swe.cotrans((lon, lat, 1.0), -obl)
            le_eq = ephem.cotrans((lon, lat, 1.0), -obl)
            diff_a = abs(se_eq[0] - le_eq[0])
            if diff_a > 180:
                diff_a = 360 - diff_a
            diff_b = abs(se_eq[1] - le_eq[1])
            if diff_a < 1e-10 and diff_b < 1e-10:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL P3 ({lon},{lat}) obl={obl}")
        except Exception as e:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: cotrans_sp comparison
print("\n=== P4: cotrans_sp SE vs LE ===")
for lon, lat in test_coords:
    for slon, slat in [(1.0, 0.0), (0.5, 0.1), (-0.5, -0.05)]:
        try:
            se_r = swe.cotrans_sp((lon, lat, 1.0, slon, slat, 0.0), -23.4393)
            le_r = ephem.cotrans_sp((lon, lat, 1.0, slon, slat, 0.0), -23.4393)
            ok = True
            for idx in range(4):
                d = abs(se_r[idx] - le_r[idx])
                if idx == 0 and d > 180:
                    d = 360 - d
                if d > 1e-8:
                    ok = False
            if ok:
                passed += 1
            else:
                failed += 1
        except Exception as e:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: J2000 equatorial
print("\n=== P5: J2000 equatorial ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.calc_ut(jd, body, 256 | 32 | 2048)
            le = ephem.calc_ut(jd, body, 256 | 32 | 2048)
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

# P6: J2000 ecliptic
print("\n=== P6: J2000 ecliptic ===")
for jd in DATES:
    for body in BODIES:
        try:
            se = swe.calc_ut(jd, body, 256 | 32)
            le = ephem.calc_ut(jd, body, 256 | 32)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 1.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# P7: Equatorial positions at different epochs
print("\n=== P7: Equatorial at various epochs ===")
for jd in [2415020.0, 2433282.0, 2451545.0, 2460000.0, 2469807.0]:
    for body in [0, 1, 4, 5]:
        try:
            se = swe.calc_ut(jd, body, 256 | 2048)
            le = ephem.calc_ut(jd, body, 256 | 2048)
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

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 72 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
