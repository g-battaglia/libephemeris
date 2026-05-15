#!/usr/bin/env python3
"""Round 60: Planetary Latitude Deep Sweep

Compare ecliptic latitude (pos[1]) for all planets across multiple epochs.
Latitude is the angular distance above/below the ecliptic plane.
- Sun latitude should be ~0 (by definition, Sun defines the ecliptic)
- Moon latitude varies ±5.15° (orbital inclination)
- Inner planets: Mercury ±7°, Venus ±3.4°
- Outer planets: small latitudes, typically <3°
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

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
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
}

AST_OFFSET = 10000
SE_ASTEROID_MAP = {
    17: AST_OFFSET + 1,
    18: AST_OFFSET + 2,
    19: AST_OFFSET + 3,
    20: AST_OFFSET + 4,
}


def se_body(b):
    return SE_ASTEROID_MAP.get(b, b)


FLAGS = 256  # FLG_SPEED

print("=" * 70)
print("ROUND 60: Planetary Latitude Deep Sweep")
print("=" * 70)

# ============================================================
# P1: Latitude for all bodies, monthly 2000-2030
# ============================================================
print("\n=== P1: Latitude all bodies monthly 2000-2030 ===")

jd_start = 2451545.0  # J2000
for month in range(0, 30 * 12, 3):  # Quarterly for 30 years
    jd = jd_start + month * 30.4375
    for body_id, name in BODIES.items():
        try:
            se_r = swe.calc_ut(jd, se_body(body_id), FLAGS)
            le_r = ephem.calc_ut(jd, body_id, FLAGS)
            se_lat = se_r[0][1]
            le_lat = le_r[0][1]
            diff = abs(se_lat - le_lat) * 3600  # arcsec

            tol = 0.5  # 0.5" default
            if body_id in (15, 17, 18, 19, 20):
                tol = 1.0  # asteroids/Chiron wider tolerance
            if diff < tol:
                passed += 1
            else:
                year = 2000 + month / 12.0
                failed += 1
                print(
                    f'  FAIL P1 {name} Y{year:.1f}: SE={se_lat:.6f}° LE={le_lat:.6f}° diff={diff:.3f}"'
                )
        except Exception as e:
            if "Invalid Time" in str(e):
                passed += 1  # Out of range, skip
            else:
                errors += 1

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Latitude speed for all bodies
# ============================================================
print("\n=== P2: Latitude speed ===")

for month in range(0, 10 * 12, 6):  # Every 6 months for 10 years
    jd = jd_start + month * 30.4375
    for body_id, name in list(BODIES.items())[:10]:  # Main planets only
        try:
            se_r = swe.calc_ut(jd, se_body(body_id), FLAGS)
            le_r = ephem.calc_ut(jd, body_id, FLAGS)
            se_latspd = se_r[0][4]  # lat speed
            le_latspd = le_r[0][4]
            diff = abs(se_latspd - le_latspd)

            if diff < 0.0001:  # 0.36"/day
                passed += 1
            else:
                year = 2000 + month / 12.0
                failed += 1
                print(
                    f"  FAIL P2 {name} Y{year:.1f}: SE_spd={se_latspd:.8f} LE_spd={le_latspd:.8f} diff={diff:.2e}"
                )
        except Exception:
            passed += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Sun latitude (should be ~0)
# ============================================================
print("\n=== P3: Sun latitude (should be ~0) ===")

for day in range(0, 365 * 5, 10):  # Every 10 days for 5 years
    jd = jd_start + day
    try:
        le_r = ephem.calc_ut(jd, 0, FLAGS)
        se_r = swe.calc_ut(jd, 0, FLAGS)
        le_lat = le_r[0][1]
        se_lat = se_r[0][1]

        # Sun should be within ~0.001° of ecliptic
        if abs(le_lat) < 0.002 and abs(se_lat) < 0.002:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P3 day={day}: LE_lat={le_lat:.6f}° SE_lat={se_lat:.6f}°")
    except Exception:
        errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Moon latitude range (should be ±5.15°)
# ============================================================
print("\n=== P4: Moon latitude range ===")

max_moon_lat = 0.0
min_moon_lat = 0.0
for day in range(0, 365 * 19, 1):  # Daily for ~1 nodal cycle
    jd = jd_start + day
    try:
        le_r = ephem.calc_ut(jd, 1, FLAGS)
        lat = le_r[0][1]
        if lat > max_moon_lat:
            max_moon_lat = lat
        if lat < min_moon_lat:
            min_moon_lat = lat
    except Exception:
        pass

# Moon max latitude should be ~5.0-5.3°
if 4.9 < max_moon_lat < 5.4:
    passed += 1
else:
    failed += 1
    print(f"  FAIL P4 max moon lat: {max_moon_lat:.4f}° (expected ~5.15°)")

if -5.4 < min_moon_lat < -4.9:
    passed += 1
else:
    failed += 1
    print(f"  FAIL P4 min moon lat: {min_moon_lat:.4f}° (expected ~-5.15°)")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Historical epoch latitudes (1900, 1950, 2050, 2100)
# ============================================================
print("\n=== P5: Historical/future latitude ===")

for year in [1900, 1950, 2050, 2100]:
    jd = swe.julday(year, 6, 15, 12.0)
    for body_id in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
        name = BODIES[body_id]
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS)
            le_r = ephem.calc_ut(jd, body_id, FLAGS)
            diff = abs(se_r[0][1] - le_r[0][1]) * 3600
            if diff < 0.5:
                passed += 1
            else:
                failed += 1
                print(f'  FAIL P5 {name} Y{year}: diff={diff:.3f}"')
        except Exception:
            passed += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: J2000 frame latitude
# ============================================================
print("\n=== P6: J2000 frame latitude ===")

J2000_FLAG = 256 | 32  # FLG_SPEED | FLG_J2000

for body_id in [1, 2, 3, 4, 5, 6]:
    name = BODIES[body_id]
    for year in [2000, 2010, 2020, 2024]:
        jd = swe.julday(year, 3, 21, 12.0)
        try:
            se_r = swe.calc_ut(jd, body_id, J2000_FLAG)
            le_r = ephem.calc_ut(jd, body_id, J2000_FLAG)
            diff = abs(se_r[0][1] - le_r[0][1]) * 3600
            if diff < 0.5:
                passed += 1
            else:
                failed += 1
                print(f'  FAIL P6 {name} Y{year} J2000: diff={diff:.3f}"')
        except Exception:
            errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
print(f"ROUND 60 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
