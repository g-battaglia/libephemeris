#!/usr/bin/env python3
"""Round 63: Barycentric Positions Deep Sweep

Compare barycentric (solar system barycenter) positions (FLG_BARYCTR).
Tests lon, lat, distance for all planets. Also BARYCTR+J2000, BARYCTR+EQUATORIAL.
"""

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
BODIES = {
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
FLG_SPEED = 256
FLG_BARYCTR = 4
FLG_J2000 = 32
FLG_EQUATORIAL = 2048

print("=" * 70)
print("ROUND 63: Barycentric Positions Deep Sweep")
print("=" * 70)

# P1: Barycentric ecliptic, monthly 2000-2025
print("\n=== P1: Barycentric ecliptic positions ===")
FLAGS = FLG_SPEED | FLG_BARYCTR
jd0 = 2451545.0
for m in range(0, 25 * 12, 2):
    jd = jd0 + m * 30.4375
    for bid, name in BODIES.items():
        try:
            se = swe.calc_ut(jd, bid, FLAGS)[0]
            le = ephem.calc_ut(jd, bid, FLAGS)[0]
            lon_d = abs(se[0] - le[0]) * 3600
            if lon_d > 180 * 3600:
                lon_d = 360 * 3600 - lon_d
            lat_d = abs(se[1] - le[1]) * 3600
            tol = 1.0 if bid in (8, 9) else 0.5  # wider for Neptune/Pluto
            if lon_d < tol and lat_d < tol:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(f'  FAIL P1 {name} Y{yr:.1f}: lon={lon_d:.3f}" lat={lat_d:.3f}"')
        except Exception as e:
            errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Barycentric speed
print("\n=== P2: Barycentric speed ===")
for m in range(0, 10 * 12, 3):
    jd = jd0 + m * 30.4375
    for bid, name in BODIES.items():
        try:
            se = swe.calc_ut(jd, bid, FLAGS)[0]
            le = ephem.calc_ut(jd, bid, FLAGS)[0]
            lon_spd = abs(se[3] - le[3])
            lat_spd = abs(se[4] - le[4])
            if lon_spd < 0.001 and lat_spd < 0.001:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(
                    f"  FAIL P2 {name} Y{yr:.1f}: lon_spd={lon_spd:.6f} lat_spd={lat_spd:.6f}"
                )
        except:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Barycentric + J2000
print("\n=== P3: Barycentric + J2000 ===")
FLAGS_J = FLG_SPEED | FLG_BARYCTR | FLG_J2000
for m in range(0, 25 * 12, 6):
    jd = jd0 + m * 30.4375
    for bid, name in BODIES.items():
        try:
            se = swe.calc_ut(jd, bid, FLAGS_J)[0]
            le = ephem.calc_ut(jd, bid, FLAGS_J)[0]
            lon_d = abs(se[0] - le[0]) * 3600
            if lon_d > 180 * 3600:
                lon_d = 360 * 3600 - lon_d
            lat_d = abs(se[1] - le[1]) * 3600
            tol = 1.0 if bid in (8, 9) else 0.5
            if lon_d < tol and lat_d < tol:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(f'  FAIL P3 {name} Y{yr:.1f}: lon={lon_d:.3f}" lat={lat_d:.3f}"')
        except:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Barycentric + Equatorial
print("\n=== P4: Barycentric + Equatorial ===")
FLAGS_EQ = FLG_SPEED | FLG_BARYCTR | FLG_EQUATORIAL
for m in range(0, 25 * 12, 6):
    jd = jd0 + m * 30.4375
    for bid, name in BODIES.items():
        try:
            se = swe.calc_ut(jd, bid, FLAGS_EQ)[0]
            le = ephem.calc_ut(jd, bid, FLAGS_EQ)[0]
            ra_d = abs(se[0] - le[0]) * 3600
            if ra_d > 180 * 3600:
                ra_d = 360 * 3600 - ra_d
            dec_d = abs(se[1] - le[1]) * 3600
            tol = 1.0 if bid in (8, 9) else 0.5
            if ra_d < tol and dec_d < tol:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(f'  FAIL P4 {name} Y{yr:.1f}: RA={ra_d:.3f}" Dec={dec_d:.3f}"')
        except:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Sun barycentric (should show barycenter-Sun offset)
print("\n=== P5: Sun barycentric position ===")
for yr in range(2000, 2026):
    jd = swe.julday(yr, 6, 15, 12.0)
    try:
        se = swe.calc_ut(jd, 0, FLAGS)[0]
        le = ephem.calc_ut(jd, 0, FLAGS)[0]
        lon_d = abs(se[0] - le[0]) * 3600
        if lon_d > 180 * 3600:
            lon_d = 360 * 3600 - lon_d
        lat_d = abs(se[1] - le[1]) * 3600
        dist_d = abs(se[2] - le[2])
        if lon_d < 0.5 and lat_d < 0.5 and dist_d < 1e-6:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P5 Sun Y{yr}: lon={lon_d:.3f}" lat={lat_d:.3f}" dist_d={dist_d:.2e}'
            )
    except:
        errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Historical barycentric (1900-2100)
print("\n=== P6: Historical barycentric ===")
for yr in range(1900, 2101, 10):
    jd = swe.julday(yr, 1, 1, 12.0)
    for bid in [0, 4, 5, 6, 7, 8, 9]:
        name = BODIES[bid]
        try:
            se = swe.calc_ut(jd, bid, FLAGS)[0]
            le = ephem.calc_ut(jd, bid, FLAGS)[0]
            lon_d = abs(se[0] - le[0]) * 3600
            if lon_d > 180 * 3600:
                lon_d = 360 * 3600 - lon_d
            tol = 2.0 if bid == 9 else 1.0 if bid == 8 else 0.5
            if lon_d < tol:
                passed += 1
            else:
                failed += 1
                print(f'  FAIL P6 {name} Y{yr}: lon={lon_d:.3f}"')
        except:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# P7: Barycentric distance
print("\n=== P7: Barycentric distance ===")
for m in range(0, 10 * 12, 6):
    jd = jd0 + m * 30.4375
    for bid in [0, 1, 4, 5, 6]:
        name = BODIES[bid]
        try:
            se = swe.calc_ut(jd, bid, FLAGS)[0]
            le = ephem.calc_ut(jd, bid, FLAGS)[0]
            dist_d = abs(se[2] - le[2])
            tol = 1e-6
            if dist_d < tol:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(f"  FAIL P7 {name} Y{yr:.1f}: dist_diff={dist_d:.2e} AU")
        except:
            errors += 1
print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
print(f"ROUND 63 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
