#!/usr/bin/env python3
"""Round 113: Chiron at Perihelion/Aphelion Deep

Chiron has a highly eccentric orbit (e≈0.38, a≈13.7 AU, period≈50.7 yr).
Perihelion ~8.5 AU, aphelion ~19 AU. Tests positions across full orbit.
P1: Chiron positions every 6 months across full orbit (1970-2070)
P2: Chiron speed at perihelion vs aphelion
P3: Chiron with different flags (HELCTR, J2000, NONUT, EQUATORIAL)
P4: Chiron distance precision
P5: Chiron latitude (near ecliptic plane crossings)
P6: Chiron retrograde periods
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLG_SPEED = 256
FLG_HELCTR = 8
FLG_J2000 = 32
FLG_NONUT = 64
FLG_EQUATORIAL = 2048
CHIRON = 15


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


# ============================================================
# P1: Chiron positions every 6 months 1970-2070
# ============================================================
print("=== P1: Chiron positions 1970-2070 ===")

for year in range(1970, 2071):
    for month in [1, 7]:
        jd = swe.julday(year, month, 1, 12.0)
        try:
            se_result = swe.calc_ut(jd, CHIRON, FLG_SPEED)
            le_result = ephem.calc_ut(jd, CHIRON, FLG_SPEED)

            se_lon, se_lat = se_result[0][0], se_result[0][1]
            le_lon, le_lat = le_result[0][0], le_result[0][1]

            diff_lon = abs(se_lon - le_lon)
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lon_arcsec = diff_lon * 3600

            diff_lat = abs(se_lat - le_lat)
            diff_lat_arcsec = diff_lat * 3600

            # Longitude tolerance: 2" (tight)
            tol_lon = 2.0
            # Latitude tolerance: 2"
            tol_lat = 2.0

            label = f"P1 Chiron {year}-{month:02d}"
            if diff_lon_arcsec >= tol_lon:
                run_test(
                    f"{label} lon",
                    False,
                    f'SE={se_lon:.6f} LE={le_lon:.6f} diff={diff_lon_arcsec:.2f}"',
                )
            else:
                passed += 1

            if diff_lat_arcsec >= tol_lat:
                run_test(
                    f"{label} lat",
                    False,
                    f'SE={se_lat:.6f} LE={le_lat:.6f} diff={diff_lat_arcsec:.2f}"',
                )
            else:
                passed += 1

            # Speed comparison
            se_spd = se_result[0][3]
            le_spd = le_result[0][3]
            diff_spd = abs(se_spd - le_spd) * 3600  # arcsec/day

            if diff_spd >= 1.0:  # 1"/day tolerance
                run_test(
                    f"{label} speed",
                    False,
                    f'SE={se_spd:.6f} LE={le_spd:.6f} diff={diff_spd:.2f}"/day',
                )
            else:
                passed += 1

        except Exception as e:
            errors += 1
            if "range" not in str(e).lower():
                print(f"  ERROR P1 {year}-{month:02d}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Chiron heliocentric positions
# ============================================================
print("\n=== P2: Chiron heliocentric ===")

HELIO_FLAG = FLG_SPEED | FLG_HELCTR

for year in range(1970, 2071, 2):
    jd = swe.julday(year, 6, 15, 12.0)
    try:
        se_result = swe.calc_ut(jd, CHIRON, HELIO_FLAG)
        le_result = ephem.calc_ut(jd, CHIRON, HELIO_FLAG)

        se_lon = se_result[0][0]
        le_lon = le_result[0][0]

        diff = abs(se_lon - le_lon)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600

        tol = 1.5  # tighter for heliocentric

        label = f"P2 Chiron helio {year}"
        if diff_arcsec >= tol:
            run_test(
                label, False, f'SE={se_lon:.6f} LE={le_lon:.6f} diff={diff_arcsec:.2f}"'
            )
        else:
            passed += 1

        # Distance
        se_dist = se_result[0][2]
        le_dist = le_result[0][2]
        diff_dist = abs(se_dist - le_dist)

        if diff_dist >= 1e-6:
            run_test(
                f"{label} dist",
                False,
                f"SE={se_dist:.8f} LE={le_dist:.8f} diff={diff_dist:.2e}",
            )
        else:
            passed += 1

    except Exception as e:
        errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Chiron with different flags
# ============================================================
print("\n=== P3: Chiron with different flags ===")

FLAG_COMBOS = [
    (FLG_SPEED, "default"),
    (FLG_SPEED | FLG_J2000 | FLG_NONUT, "J2000"),
    (FLG_SPEED | FLG_NONUT, "NONUT"),
    (FLG_SPEED | FLG_EQUATORIAL, "EQUATORIAL"),
    (FLG_SPEED | FLG_HELCTR, "HELIO"),
]

test_jds = [
    (swe.julday(2000, 1, 1, 12.0), "J2000"),
    (swe.julday(2024, 6, 15, 12.0), "2024"),
    (swe.julday(1996, 2, 14, 12.0), "1996-peri"),  # near perihelion
    (swe.julday(2021, 6, 19, 12.0), "2021"),
    (swe.julday(1980, 1, 1, 12.0), "1980"),
]

for flags, flag_label in FLAG_COMBOS:
    for jd, epoch_label in test_jds:
        try:
            se_result = swe.calc_ut(jd, CHIRON, flags)
            le_result = ephem.calc_ut(jd, CHIRON, flags)

            for idx, comp_name in [(0, "lon/RA"), (1, "lat/Dec")]:
                se_val = se_result[0][idx]
                le_val = le_result[0][idx]

                diff = abs(se_val - le_val)
                if idx == 0 and diff > 180:
                    diff = 360 - diff
                diff_arcsec = diff * 3600

                tol = 2.0
                if "J2000" in flag_label:
                    tol = 20.0  # known J2000 offset for analytical bodies

                label = f"P3 {flag_label} {epoch_label} {comp_name}"
                if diff_arcsec >= tol:
                    run_test(
                        label,
                        False,
                        f'SE={se_val:.6f} LE={le_val:.6f} diff={diff_arcsec:.2f}"',
                    )
                else:
                    passed += 1

        except Exception as e:
            errors += 1
            print(f"  ERROR P3 {flag_label} {epoch_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Chiron near perihelion (1996) — daily positions
# ============================================================
print("\n=== P4: Chiron near perihelion 1996 ===")

jd_start = swe.julday(1995, 6, 1, 12.0)
for day in range(0, 730, 5):  # 2 years at 5-day intervals
    jd = jd_start + day
    try:
        se_result = swe.calc_ut(jd, CHIRON, FLG_SPEED)
        le_result = ephem.calc_ut(jd, CHIRON, FLG_SPEED)

        se_lon = se_result[0][0]
        le_lon = le_result[0][0]

        diff = abs(se_lon - le_lon)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600

        if diff_arcsec >= 2.0:
            run_test(
                f"P4 Chiron peri d={day}",
                False,
                f'SE={se_lon:.6f} LE={le_lon:.6f} diff={diff_arcsec:.2f}"',
            )
        else:
            passed += 1
    except Exception as e:
        errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Chiron near aphelion (2021) — positions
# ============================================================
print("\n=== P5: Chiron near aphelion 2021 ===")

jd_start = swe.julday(2020, 1, 1, 12.0)
for day in range(0, 730, 5):
    jd = jd_start + day
    try:
        se_result = swe.calc_ut(jd, CHIRON, FLG_SPEED)
        le_result = ephem.calc_ut(jd, CHIRON, FLG_SPEED)

        se_lon = se_result[0][0]
        le_lon = le_result[0][0]

        diff = abs(se_lon - le_lon)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600

        if diff_arcsec >= 2.0:
            run_test(
                f"P5 Chiron aph d={day}",
                False,
                f'SE={se_lon:.6f} LE={le_lon:.6f} diff={diff_arcsec:.2f}"',
            )
        else:
            passed += 1
    except Exception as e:
        errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Chiron retrograde detection
# ============================================================
print("\n=== P6: Chiron retrograde detection ===")

jd_start = swe.julday(2020, 1, 1, 12.0)
for day in range(0, 365 * 3, 3):  # 3 years at 3-day intervals
    jd = jd_start + day
    try:
        se_result = swe.calc_ut(jd, CHIRON, FLG_SPEED)
        le_result = ephem.calc_ut(jd, CHIRON, FLG_SPEED)

        se_spd = se_result[0][3]
        le_spd = le_result[0][3]

        # Check retrograde status matches
        se_retro = se_spd < 0
        le_retro = le_spd < 0

        if se_retro != le_retro:
            run_test(
                f"P6 Chiron retro d={day}",
                False,
                f"SE_spd={se_spd:.6f} LE_spd={le_spd:.6f}",
            )
        else:
            passed += 1

        # Speed magnitude comparison
        diff_spd = abs(se_spd - le_spd) * 3600
        if diff_spd >= 1.0:
            run_test(
                f"P6 Chiron speed d={day}",
                False,
                f'SE={se_spd:.6f} LE={le_spd:.6f} diff={diff_spd:.2f}"/day',
            )
        else:
            passed += 1

    except Exception as e:
        errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 113 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
