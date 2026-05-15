#!/usr/bin/env python3
"""Round 71: Topocentric Moon Parallax"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256 | 32768  # FLG_SPEED | FLG_TOPOCTR

LOCATIONS = [
    (12.50, 41.90, 50.0, "Rome"),
    (-73.97, 40.78, 10.0, "NYC"),
    (139.69, 35.69, 40.0, "Tokyo"),
    (151.21, -33.87, 20.0, "Sydney"),
    (0.0, 0.0, 0.0, "Null"),
    (18.07, 59.33, 15.0, "Stockholm"),
    (-43.17, -22.91, 11.0, "Rio"),
    (77.21, 28.61, 216.0, "Delhi"),
]

DATES = [2451545.0 + i * 30 for i in range(40)]  # 40 months

print("=" * 70)
print("ROUND 71: Topocentric Positions (Moon Parallax Focus)")
print("=" * 70)

# P1: Topocentric Moon positions
print("\n=== P1: Topocentric Moon ===")
for lon, lat, alt, name in LOCATIONS:
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)
    loc_pass = loc_fail = 0
    for jd in DATES:
        try:
            se = swe.calc_ut(jd, 1, FLAGS)
            le = ephem.calc_ut(jd, 1, FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            diff_as = diff * 3600
            if diff_as < 1.0:
                passed += 1
                loc_pass += 1
            else:
                failed += 1
                loc_fail += 1
                if loc_fail <= 2:
                    print(
                        f'  FAIL {name} jd={jd:.0f}: SE={se[0][0]:.6f} LE={le[0][0]:.6f} d={diff_as:.2f}"'
                    )
        except Exception as e:
            errors += 1
    if loc_fail > 0:
        print(f"  {name}: {loc_pass}/{loc_pass + loc_fail}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Topocentric Sun
print("\n=== P2: Topocentric Sun ===")
for lon, lat, alt, name in LOCATIONS:
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)
    for jd in DATES:
        try:
            se = swe.calc_ut(jd, 0, FLAGS)
            le = ephem.calc_ut(jd, 0, FLAGS)
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

# P3: Topocentric planets
print("\n=== P3: Topocentric planets ===")
for lon, lat, alt, name in LOCATIONS[:4]:
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)
    for jd in DATES[:10]:
        for body in [2, 3, 4, 5, 6]:
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

# P4: Moon parallax magnitude check (should be ~0.5-1°)
print("\n=== P4: Moon parallax magnitude ===")
swe.set_topo(12.50, 41.90, 50.0)
ephem.set_topo(12.50, 41.90, 50.0)
for jd in DATES[:20]:
    try:
        geo = ephem.calc_ut(jd, 1, 256)  # geocentric
        topo = ephem.calc_ut(jd, 1, FLAGS)  # topocentric
        parallax = abs(geo[0][0] - topo[0][0])
        if parallax > 180:
            parallax = 360 - parallax
        # Moon parallax should be between 0.1° and 1.5°
        if 0.0 < parallax < 1.5:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P4 jd={jd:.0f}: parallax={parallax:.4f}°")
    except:
        errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Topocentric lat comparison
print("\n=== P5: Topocentric latitude ===")
for lon, lat, alt, name in LOCATIONS[:4]:
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)
    for jd in DATES[:10]:
        for body in [0, 1]:
            try:
                se = swe.calc_ut(jd, body, FLAGS)
                le = ephem.calc_ut(jd, body, FLAGS)
                diff = abs(se[0][1] - le[0][1]) * 3600
                if diff < 1.0:
                    passed += 1
                else:
                    failed += 1
            except:
                errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 71 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
