#!/usr/bin/env python3
"""Round 65: Gauquelin Sectors with Different Methods

Compare gauquelin_sector() across methods 0-5, multiple planets,
locations, and dates. Methods 0-1 use hour-angle approximation,
methods 2-5 use actual rise/set times.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0

# Bodies to test
BODIES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  # Sun through Pluto
BODY_NAMES = [
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

# Locations: (lon, lat, alt, name)
LOCATIONS = [
    (2.35, 48.85, 0.0, "Paris"),
    (-73.97, 40.78, 0.0, "New York"),
    (139.69, 35.69, 0.0, "Tokyo"),
    (12.50, 41.90, 0.0, "Rome"),
    (-43.17, -22.91, 0.0, "Rio"),
    (151.21, -33.87, 0.0, "Sydney"),
    (0.0, 0.0, 0.0, "Null Island"),
    (18.07, 59.33, 0.0, "Stockholm"),
    (24.94, 60.17, 0.0, "Helsinki"),
]

# Test dates
DATES = [
    2451545.0,  # J2000 - 2000-01-01 12:00 UT
    2451635.0,  # 2000-04-01
    2451727.0,  # 2000-07-02
    2451818.0,  # 2000-10-01
    2460000.0,  # 2023-02-25
    2460400.0,  # 2024-03-31
]

print("=" * 70)
print("ROUND 65: Gauquelin Sectors with Different Methods")
print("=" * 70)

# ============================================================
# P1: Method 0 (with latitude) — all bodies, locations, dates
# ============================================================
print("\n=== P1: Method 0 (with latitude) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS:
        geopos = (lon, lat, alt)
        for body in BODIES:
            try:
                se_result = swe.gauquelin_sector(jd, body, 0, geopos, 1013.25, 15.0)
                le_result = ephem.gauquelin_sector(
                    jd, body, 0, geopos, 1013.25, 15.0
                )
                diff = abs(se_result - le_result)
                # Handle wrap (sector 36 → 1)
                if diff > 18:
                    diff = 36 - diff
                if diff < 0.15:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P1 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
                if "circumpolar" not in str(e).lower():
                    print(f"  ERR P1 {loc_name} body={BODY_NAMES[body]} jd={jd}: {e}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Method 1 (without latitude) — all bodies, locations, dates
# ============================================================
print("\n=== P2: Method 1 (without latitude) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS:
        geopos = (lon, lat, alt)
        for body in BODIES:
            try:
                se_result = swe.gauquelin_sector(jd, body, 1, geopos, 1013.25, 15.0)
                le_result = ephem.gauquelin_sector(
                    jd, body, 1, geopos, 1013.25, 15.0
                )
                diff = abs(se_result - le_result)
                if diff > 18:
                    diff = 36 - diff
                if diff < 0.15:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P2 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
                if "circumpolar" not in str(e).lower():
                    print(f"  ERR P2 {loc_name} body={BODY_NAMES[body]} jd={jd}: {e}")
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Method 2 (disc center, no refraction) — Sun, Moon, Mars
# ============================================================
print("\n=== P3: Method 2 (disc center, no refraction) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS:
        if abs(lat) > 60:
            continue  # Skip high latitudes for rise/set methods
        geopos = (lon, lat, alt)
        for body in [0, 1, 4]:  # Sun, Moon, Mars
            try:
                se_result = swe.gauquelin_sector(jd, body, 2, geopos, 1013.25, 15.0)
                le_result = ephem.gauquelin_sector(
                    jd, body, 2, geopos, 1013.25, 15.0
                )
                diff = abs(se_result - le_result)
                if diff > 18:
                    diff = 36 - diff
                if diff < 0.5:  # Rise/set methods have more inherent variation
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P3 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
                if "circumpolar" not in str(e).lower():
                    print(f"  ERR P3 {loc_name} body={BODY_NAMES[body]} jd={jd}: {e}")
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Method 3 (disc center, with refraction) — Sun, Moon, Mars
# ============================================================
print("\n=== P4: Method 3 (disc center, with refraction) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS:
        if abs(lat) > 60:
            continue
        geopos = (lon, lat, alt)
        for body in [0, 1, 4]:
            try:
                se_result = swe.gauquelin_sector(jd, body, 3, geopos, 1013.25, 15.0)
                le_result = ephem.gauquelin_sector(
                    jd, body, 3, geopos, 1013.25, 15.0
                )
                diff = abs(se_result - le_result)
                if diff > 18:
                    diff = 36 - diff
                if diff < 0.5:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P4 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
                if "circumpolar" not in str(e).lower():
                    print(f"  ERR P4 {loc_name} body={BODY_NAMES[body]} jd={jd}: {e}")
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Method 4 (disc edge, no refraction) — Sun, Moon
# ============================================================
print("\n=== P5: Method 4 (disc edge, no refraction) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS:
        if abs(lat) > 55:
            continue
        geopos = (lon, lat, alt)
        for body in [0, 1]:
            try:
                se_result = swe.gauquelin_sector(jd, body, 4, geopos, 1013.25, 15.0)
                le_result = ephem.gauquelin_sector(
                    jd, body, 4, geopos, 1013.25, 15.0
                )
                diff = abs(se_result - le_result)
                if diff > 18:
                    diff = 36 - diff
                if diff < 0.5:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P5 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
                if "circumpolar" not in str(e).lower():
                    print(f"  ERR P5 {loc_name} body={BODY_NAMES[body]} jd={jd}: {e}")
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: Method 5 (disc edge, with refraction) — Sun, Moon
# ============================================================
print("\n=== P6: Method 5 (disc edge, with refraction) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS:
        if abs(lat) > 55:
            continue
        geopos = (lon, lat, alt)
        for body in [0, 1]:
            try:
                se_result = swe.gauquelin_sector(jd, body, 5, geopos, 1013.25, 15.0)
                le_result = ephem.gauquelin_sector(
                    jd, body, 5, geopos, 1013.25, 15.0
                )
                diff = abs(se_result - le_result)
                if diff > 18:
                    diff = 36 - diff
                if diff < 0.5:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P6 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
                    )
            except Exception as e:
                errors += 1
                if "circumpolar" not in str(e).lower():
                    print(f"  ERR P6 {loc_name} body={BODY_NAMES[body]} jd={jd}: {e}")
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Sector range validation — all methods
# ============================================================
print("\n=== P7: Sector range validation ===")
for method in range(6):
    for jd in [2451545.0, 2460000.0]:
        for lon, lat, alt, loc_name in LOCATIONS:
            if abs(lat) > 55 and method >= 2:
                continue
            geopos = (lon, lat, alt)
            for body in [0, 1, 2, 4]:
                try:
                    result = ephem.gauquelin_sector(
                        jd, body, method, geopos, 1013.25, 15.0
                    )
                    if 1.0 <= result < 37.0:
                        passed += 1
                    else:
                        failed += 1
                        print(
                            f"  FAIL P7 method={method} {loc_name} body={BODY_NAMES[body]}: "
                            f"sector={result:.4f} (out of range [1, 37))"
                        )
                except Exception as e:
                    errors += 1
print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P8: Method 0 vs Method 1 consistency (with/without latitude)
# ============================================================
print("\n=== P8: Method 0 vs 1 consistency (lat matters for Moon) ===")
for jd in DATES:
    for lon, lat, alt, loc_name in LOCATIONS[:5]:
        geopos = (lon, lat, alt)
        for body in BODIES:
            try:
                le_m0 = ephem.gauquelin_sector(jd, body, 0, geopos, 1013.25, 15.0)
                le_m1 = ephem.gauquelin_sector(jd, body, 1, geopos, 1013.25, 15.0)
                se_m0 = swe.gauquelin_sector(jd, body, 0, geopos, 1013.25, 15.0)
                se_m1 = swe.gauquelin_sector(jd, body, 1, geopos, 1013.25, 15.0)
                # The differences between method 0/1 should be the same in both
                le_diff = le_m0 - le_m1
                se_diff = se_m0 - se_m1
                delta = abs(le_diff - se_diff)
                if delta < 0.15:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P8 {loc_name} body={BODY_NAMES[body]} jd={jd}: "
                        f"LE_diff={le_diff:.4f} SE_diff={se_diff:.4f} delta={delta:.4f}"
                    )
            except Exception as e:
                errors += 1
print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P9: Dense time sweep — sector continuity for Sun
# ============================================================
print("\n=== P9: Dense time sweep for Sun sectors ===")
jd_base = 2451545.0
geopos = (12.50, 41.90, 0.0)  # Rome
prev_le = None
prev_se = None
for i in range(48):  # 48 half-hours = 1 full day
    jd = jd_base + i / 48.0
    try:
        se_result = swe.gauquelin_sector(jd, 0, 0, geopos, 1013.25, 15.0)
        le_result = ephem.gauquelin_sector(jd, 0, 0, geopos, 1013.25, 15.0)
        diff = abs(se_result - le_result)
        if diff > 18:
            diff = 36 - diff
        if diff < 0.15:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P9 i={i} jd={jd:.4f}: SE={se_result:.4f} LE={le_result:.4f} diff={diff:.4f}"
            )
    except Exception as e:
        errors += 1
print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 65 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
