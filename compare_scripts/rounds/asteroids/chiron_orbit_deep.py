#!/usr/bin/env python3
"""Round 84: Chiron Full Orbit Deep Verification

Chiron has a ~50.7 year orbit with high eccentricity (0.38) and inclination (6.9°).
Tests positions across multiple complete orbits, including perihelion/aphelion,
node crossings, and with various flag combinations.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLG_SWIEPH = 2
FLG_SPEED = 256
FLG_HELCTR = 8
FLG_J2000 = 32
FLG_NONUT = 64
FLG_EQUATORIAL = 2048
CHIRON = 15

print("=" * 70)
print("ROUND 84: Chiron Full Orbit Deep Verification")
print("=" * 70)


def check(label, se_val, le_val, tol_arcsec):
    global passed, failed
    diff = abs(se_val - le_val) * 3600.0
    if diff > 180 * 3600:
        diff = 360 * 3600 - diff
    if diff < tol_arcsec:
        passed += 1
    else:
        failed += 1
        print(f'  FAIL {label}: SE={se_val:.6f} LE={le_val:.6f} diff={diff:.2f}"')


# ============================================================
# P1: Geocentric ecliptic positions across 100 years (monthly)
# ============================================================
print("\n=== P1: Geocentric ecliptic longitude (monthly 1950-2050) ===")

flags = FLG_SWIEPH | FLG_SPEED
for year in range(1950, 2051, 2):
    for month in [1, 7]:
        jd = swe.julday(year, month, 1, 12.0)
        label = f"Chiron lon {year}-{month:02d}"
        try:
            se = swe.calc_ut(jd, CHIRON, flags)
            le = ephem.calc_ut(jd, CHIRON, flags)
            check(label, se[0][0], le[0][0], 1.0)
        except Exception as e:
            errors += 1
            if errors <= 3:
                print(f"  ERROR {label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Latitude
# ============================================================
print("\n=== P2: Geocentric ecliptic latitude ===")

for year in range(1950, 2051, 2):
    for month in [1, 7]:
        jd = swe.julday(year, month, 1, 12.0)
        label = f"Chiron lat {year}-{month:02d}"
        try:
            se = swe.calc_ut(jd, CHIRON, flags)
            le = ephem.calc_ut(jd, CHIRON, flags)
            check(label, se[0][1], le[0][1], 1.0)
        except Exception as e:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Distance (AU)
# ============================================================
print("\n=== P3: Distance (AU) ===")

for year in range(1950, 2051, 5):
    jd = swe.julday(year, 1, 1, 12.0)
    label = f"Chiron dist {year}"
    try:
        se = swe.calc_ut(jd, CHIRON, flags)
        le = ephem.calc_ut(jd, CHIRON, flags)
        ratio = se[0][2] / le[0][2] if le[0][2] != 0 else 999
        if abs(ratio - 1.0) < 0.0001:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL {label}: SE={se[0][2]:.8f} LE={le[0][2]:.8f} ratio={ratio:.8f}"
            )
    except Exception as e:
        errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Speed (longitude speed)
# ============================================================
print("\n=== P4: Longitude speed ===")

for year in range(1950, 2051, 2):
    jd = swe.julday(year, 6, 15, 12.0)
    label = f"Chiron speed {year}"
    try:
        se = swe.calc_ut(jd, CHIRON, flags)
        le = ephem.calc_ut(jd, CHIRON, flags)
        diff = abs(se[0][3] - le[0][3]) * 3600.0
        if diff < 1.0:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL {label}: SE={se[0][3]:.6f} LE={le[0][3]:.6f} diff={diff:.2f}"/day'
            )
    except Exception as e:
        errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Heliocentric positions
# ============================================================
print("\n=== P5: Heliocentric positions ===")

hflags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
for year in range(1950, 2051, 5):
    jd = swe.julday(year, 1, 1, 12.0)
    label = f"Chiron helio {year}"
    try:
        se = swe.calc_ut(jd, CHIRON, hflags)
        le = ephem.calc_ut(jd, CHIRON, hflags)
        check(f"{label} lon", se[0][0], le[0][0], 1.0)
        check(f"{label} lat", se[0][1], le[0][1], 1.0)
    except Exception as e:
        errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: J2000 frame
# ============================================================
print("\n=== P6: J2000 ecliptic frame ===")

j2kflags = FLG_SWIEPH | FLG_SPEED | FLG_J2000 | FLG_NONUT
for year in range(1950, 2051, 5):
    jd = swe.julday(year, 1, 1, 12.0)
    label = f"Chiron J2000 {year}"
    try:
        se = swe.calc_ut(jd, CHIRON, j2kflags)
        le = ephem.calc_ut(jd, CHIRON, j2kflags)
        check(f"{label} lon", se[0][0], le[0][0], 1.0)
        check(f"{label} lat", se[0][1], le[0][1], 1.0)
    except Exception as e:
        errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Equatorial coordinates (RA/Dec)
# ============================================================
print("\n=== P7: Equatorial coordinates ===")

eqflags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
for year in range(1950, 2051, 5):
    jd = swe.julday(year, 1, 1, 12.0)
    label = f"Chiron eq {year}"
    try:
        se = swe.calc_ut(jd, CHIRON, eqflags)
        le = ephem.calc_ut(jd, CHIRON, eqflags)
        check(f"{label} RA", se[0][0], le[0][0], 1.5)
        check(f"{label} Dec", se[0][1], le[0][1], 1.5)
    except Exception as e:
        errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P8: Around perihelion (1996, ~8.45 AU)
# ============================================================
print("\n=== P8: Near perihelion (1996) ===")

for month in range(1, 13):
    jd = swe.julday(1996, month, 15, 12.0)
    label = f"Chiron perihelion 1996-{month:02d}"
    try:
        se = swe.calc_ut(jd, CHIRON, flags)
        le = ephem.calc_ut(jd, CHIRON, flags)
        check(f"{label} lon", se[0][0], le[0][0], 1.0)
        check(f"{label} lat", se[0][1], le[0][1], 1.0)
    except Exception as e:
        errors += 1

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P9: Around aphelion (~2021, ~18.8 AU)
# ============================================================
print("\n=== P9: Near aphelion (2021) ===")

for month in range(1, 13):
    jd = swe.julday(2021, month, 15, 12.0)
    label = f"Chiron aphelion 2021-{month:02d}"
    try:
        se = swe.calc_ut(jd, CHIRON, flags)
        le = ephem.calc_ut(jd, CHIRON, flags)
        check(f"{label} lon", se[0][0], le[0][0], 1.0)
        check(f"{label} lat", se[0][1], le[0][1], 1.0)
    except Exception as e:
        errors += 1

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P10: Retrograde periods (Chiron retrogrades ~5 months/year)
# ============================================================
print("\n=== P10: Retrograde detection ===")

retro_match = 0
retro_mismatch = 0
for year in range(1980, 2030):
    for month in range(1, 13):
        jd = swe.julday(year, month, 15, 12.0)
        try:
            se = swe.calc_ut(jd, CHIRON, flags)
            le = ephem.calc_ut(jd, CHIRON, flags)
            se_retro = se[0][3] < 0
            le_retro = le[0][3] < 0
            if se_retro == le_retro:
                retro_match += 1
            else:
                retro_mismatch += 1
                if retro_mismatch <= 3:
                    print(
                        f"  FAIL retro {year}-{month:02d}: SE_retro={se_retro} LE_retro={le_retro} SE_spd={se[0][3]:.6f} LE_spd={le[0][3]:.6f}"
                    )
        except Exception as e:
            errors += 1

passed += retro_match
failed += retro_mismatch
print(f"  Retrograde match: {retro_match}, mismatch: {retro_mismatch}")
print(f"  After P10: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# Summary
# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 84 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
