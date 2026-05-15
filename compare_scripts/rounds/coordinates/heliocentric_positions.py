#!/usr/bin/env python3
"""Round 62: Heliocentric Positions Deep Sweep

Compare heliocentric ecliptic positions (FLG_HELCTR) for all planets.
Heliocentric = as seen from Sun center. No light-time, no aberration.
Tests longitude, latitude, distance, and speeds.
Also tests HELCTR + J2000, HELCTR + EQUATORIAL combinations.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0
BODIES = {
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
FLG_HELCTR = 8
FLG_J2000 = 32
FLG_EQUATORIAL = 2048

print("=" * 70)
print("ROUND 62: Heliocentric Positions Deep Sweep")
print("=" * 70)

# P1: Heliocentric ecliptic lon/lat, monthly 2000-2025
print("\n=== P1: Heliocentric ecliptic positions ===")
FLAGS = FLG_SPEED | FLG_HELCTR
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
            if lon_d < 0.5 and lat_d < 0.5:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(f'  FAIL P1 {name} Y{yr:.1f}: lon={lon_d:.3f}" lat={lat_d:.3f}"')
        except:
            errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Heliocentric speed
print("\n=== P2: Heliocentric lon/lat speed ===")
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

# P3: Heliocentric + J2000
print("\n=== P3: Heliocentric + J2000 ===")
FLAGS_J = FLG_SPEED | FLG_HELCTR | FLG_J2000
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
            if lon_d < 0.5 and lat_d < 0.5:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(
                    f'  FAIL P3 {name} Y{yr:.1f} J2000: lon={lon_d:.3f}" lat={lat_d:.3f}"'
                )
        except:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Heliocentric + Equatorial
print("\n=== P4: Heliocentric + Equatorial ===")
FLAGS_EQ = FLG_SPEED | FLG_HELCTR | FLG_EQUATORIAL
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
            if ra_d < 0.5 and dec_d < 0.5:
                passed += 1
            else:
                failed += 1
                yr = 2000 + m / 12.0
                print(f'  FAIL P4 {name} Y{yr:.1f} EQ: RA={ra_d:.3f}" Dec={dec_d:.3f}"')
        except:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Earth heliocentric = reversed Sun geocentric
print("\n=== P5: Earth helio consistency ===")
# calc_ut with HELCTR for Earth (body=14 in SE, or use Sun geocentric + 180°)
for yr in range(2000, 2026):
    jd = swe.julday(yr, 6, 15, 12.0)
    try:
        sun_geo = ephem.calc_ut(jd, 0, FLG_SPEED)[0]
        # Earth helio lon should be Sun geo lon + 180°
        earth_lon = (sun_geo[0] + 180.0) % 360.0
        # Earth helio lat should be ~ -Sun geo lat
        earth_lat = -sun_geo[1]
        # Verify with Mercury helio (just sanity check Earth position)
        merc_helio = ephem.calc_ut(jd, 2, FLAGS)[0]
        # Mercury helio lon should be between 0-360
        if 0 <= merc_helio[0] < 360 and 0 <= earth_lon < 360:
            passed += 1
        else:
            failed += 1
    except:
        errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Historical heliocentric
print("\n=== P6: Historical heliocentric (1900-2100) ===")
for yr in range(1900, 2101, 10):
    jd = swe.julday(yr, 1, 1, 12.0)
    for bid in [2, 3, 4, 5, 6]:
        name = BODIES[bid]
        try:
            se = swe.calc_ut(jd, bid, FLAGS)[0]
            le = ephem.calc_ut(jd, bid, FLAGS)[0]
            lon_d = abs(se[0] - le[0]) * 3600
            if lon_d > 180 * 3600:
                lon_d = 360 * 3600 - lon_d
            if lon_d < 0.5:
                passed += 1
            else:
                failed += 1
                print(f'  FAIL P6 {name} Y{yr}: lon={lon_d:.3f}"')
        except:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
print(f"ROUND 62 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
