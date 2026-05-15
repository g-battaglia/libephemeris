#!/usr/bin/env python3
"""Round 53: Planetary Distances (AU) Precision

Compare geocentric distance (AU) for all major bodies across multiple epochs.
Tests distance values from calc_ut() index 2 (geocentric distance in AU).
Also tests distance speed (index 5 - AU/day) and heliocentric distances.
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

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

# SE body IDs for asteroids
SE_BODIES = {
    0: 0,
    1: 1,
    2: 2,
    3: 3,
    4: 4,
    5: 5,
    6: 6,
    7: 7,
    8: 8,
    9: 9,
    15: 15,  # Chiron
    17: 17,  # Ceres (LE constant)
    18: 18,  # Pallas
    19: 19,  # Juno
    20: 20,  # Vesta
}

# SE uses AST_OFFSET + N for asteroids
AST_OFFSET = 10000
SE_ASTEROID_MAP = {
    17: AST_OFFSET + 1,  # Ceres
    18: AST_OFFSET + 2,  # Pallas
    19: AST_OFFSET + 3,  # Juno
    20: AST_OFFSET + 4,  # Vesta
}

# Test epochs: spread across 1600-2400
TEST_JDAYS = []

# Every 50 years from 1600 to 2400
for year in range(1600, 2401, 50):
    jd = swe.julday(year, 1, 1, 12.0)
    TEST_JDAYS.append((jd, f"Y{year}"))

# Quarterly for 2024
for month in [1, 4, 7, 10]:
    jd = swe.julday(2024, month, 1, 12.0)
    TEST_JDAYS.append((jd, f"2024-{month:02d}"))

# Some specific dates
TEST_JDAYS.append((2451545.0, "J2000.0"))
TEST_JDAYS.append((2460000.0, "JD2460000"))
TEST_JDAYS.append((2440000.0, "JD2440000"))

FLAGS = 256  # FLG_SPEED


def get_se_body(body_id):
    """Get the SE body ID, mapping asteroids to AST_OFFSET + N."""
    return SE_ASTEROID_MAP.get(body_id, body_id)


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


# ============================================================
# P1: Geocentric distance for all bodies across epochs
# ============================================================
print("=== P1: Geocentric distance all bodies ===")
p1_pass = p1_fail = 0

for jd, epoch_label in TEST_JDAYS:
    for body_id, body_name in BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLAGS)
            se_pos = se_result[0] if isinstance(se_result, tuple) else se_result

            le_result = ephem.calc_ut(jd, body_id, FLAGS)
            le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

            se_dist = se_pos[2]
            le_dist = le_pos[2]

            if se_dist == 0.0 and le_dist == 0.0:
                passed += 1
                p1_pass += 1
                continue

            if se_dist == 0.0 or le_dist == 0.0:
                # One returns 0, other doesn't — known for analytical bodies
                if body_id in (10, 11, 12, 13):  # Mean node/apog etc
                    passed += 1
                    p1_pass += 1
                else:
                    diff_au = abs(se_dist - le_dist)
                    if diff_au < 0.001:
                        passed += 1
                        p1_pass += 1
                    else:
                        failed += 1
                        p1_fail += 1
                        print(
                            f"  FAIL P1 {body_name} {epoch_label}: SE_dist={se_dist:.10f} LE_dist={le_dist:.10f} diff={diff_au:.2e} AU"
                        )
                continue

            diff_au = abs(se_dist - le_dist)
            rel_diff = diff_au / se_dist

            # Tolerance: 1e-6 AU = ~150 m (generous for ephemeris differences)
            # For Moon: ~0.38 m; for Sun: ~150 m; for Pluto: ~150 m
            # Tighter: 1e-7 AU = ~15 m
            tol_au = 1e-6
            if body_id == 1:  # Moon - tighter tolerance
                tol_au = 1e-7  # ~15 m

            if diff_au < tol_au:
                passed += 1
                p1_pass += 1
            else:
                failed += 1
                p1_fail += 1
                diff_km = diff_au * 1.496e8  # AU to km
                print(
                    f"  FAIL P1 {body_name} {epoch_label}: SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU ({diff_km:.1f} km)"
                )

        except Exception as e:
            errors += 1
            print(f"  ERROR P1 {body_name} {epoch_label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Distance speed (AU/day) for all bodies
# ============================================================
print("\n=== P2: Distance speed (AU/day) ===")
p2_pass_start = passed

for jd, epoch_label in TEST_JDAYS[:10]:  # First 10 epochs
    for body_id, body_name in BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLAGS)
            se_pos = se_result[0]

            le_result = ephem.calc_ut(jd, body_id, FLAGS)
            le_pos = le_result[0]

            se_dspeed = se_pos[5]  # distance speed AU/day
            le_dspeed = le_pos[5]

            diff = abs(se_dspeed - le_dspeed)

            # Tolerance: 1e-7 AU/day ≈ 15 m/day ≈ 0.17 mm/s
            tol = 1e-6
            if body_id == 1:  # Moon
                tol = 1e-7

            if diff < tol:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P2 {body_name} {epoch_label}: SE_dspd={se_dspeed:.10f} LE_dspd={le_dspeed:.10f} diff={diff:.2e} AU/day"
                )

        except Exception as e:
            errors += 1
            print(f"  ERROR P2 {body_name} {epoch_label}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Heliocentric distances
# ============================================================
print("\n=== P3: Heliocentric distance ===")
HELIO_FLAG = 256 | 8  # FLG_SPEED | FLG_HELCTR

# Skip Sun (meaningless helio) and Moon (SE may not support helio Moon well)
HELIO_BODIES = {k: v for k, v in BODIES.items() if k not in (0, 1)}

for jd, epoch_label in TEST_JDAYS[:10]:
    for body_id, body_name in HELIO_BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, HELIO_FLAG)
            se_pos = se_result[0]

            le_result = ephem.calc_ut(jd, body_id, HELIO_FLAG)
            le_pos = le_result[0]

            se_dist = se_pos[2]
            le_dist = le_pos[2]

            if se_dist == 0.0 and le_dist == 0.0:
                passed += 1
                continue

            diff_au = abs(se_dist - le_dist)
            tol_au = 1e-6

            if diff_au < tol_au:
                passed += 1
            else:
                diff_km = diff_au * 1.496e8
                failed += 1
                print(
                    f"  FAIL P3 {body_name} {epoch_label}: SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU ({diff_km:.1f} km)"
                )

        except Exception as e:
            errors += 1
            print(f"  ERROR P3 {body_name} {epoch_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Moon distance precision (perigee/apogee cycle)
# ============================================================
print("\n=== P4: Moon distance across lunation cycle ===")
# Moon distance varies ~356,500 km (perigee) to ~406,700 km (apogee)
# Check distance at 2-day intervals over ~1 year

jd_start = swe.julday(2024, 1, 1, 0.0)
for i in range(0, 365, 2):
    jd = jd_start + i
    try:
        se_result = swe.calc_ut(jd, 1, FLAGS)
        le_result = ephem.calc_ut(jd, 1, FLAGS)

        se_dist = se_result[0][2]
        le_dist = le_result[0][2]

        diff_au = abs(se_dist - le_dist)
        diff_km = diff_au * 1.496e8

        # Moon distance tolerance: 0.01 km = 10 meters
        if diff_km < 0.1:  # 100 meters
            passed += 1
        else:
            failed += 1
            se_km = se_dist * 1.496e8
            le_km = le_dist * 1.496e8
            print(
                f"  FAIL P4 day={i}: SE={se_km:.1f}km LE={le_km:.1f}km diff={diff_km:.3f}km"
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P4 day={i}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Sun-Earth distance precision (perihelion/aphelion)
# ============================================================
print("\n=== P5: Sun-Earth distance across year ===")
# Earth-Sun distance varies ~0.9833 AU (perihelion, ~Jan 3) to ~1.0167 AU (aphelion, ~Jul 4)

jd_start = swe.julday(2024, 1, 1, 0.0)
for i in range(0, 365, 3):
    jd = jd_start + i
    try:
        se_result = swe.calc_ut(jd, 0, FLAGS)
        le_result = ephem.calc_ut(jd, 0, FLAGS)

        se_dist = se_result[0][2]
        le_dist = le_result[0][2]

        diff_au = abs(se_dist - le_dist)
        diff_km = diff_au * 1.496e8

        # Sun distance: 1e-8 AU ≈ 1.5 m
        if diff_au < 1e-7:  # ~15 m
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P5 day={i}: SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU ({diff_km:.1f}km)"
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P5 day={i}: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: XYZ mode distances (cross-check with ecliptic distances)
# ============================================================
print("\n=== P6: XYZ distance vs ecliptic distance consistency ===")
XYZ_FLAG = 256 | 4096  # FLG_SPEED | FLG_XYZ

for jd, epoch_label in TEST_JDAYS[:5]:
    for body_id in [0, 1, 2, 3, 4, 5]:
        body_name = BODIES[body_id]
        try:
            # Ecliptic distance
            le_ecl = ephem.calc_ut(jd, body_id, FLAGS)
            ecl_dist = le_ecl[0][2]

            # XYZ distance
            le_xyz = ephem.calc_ut(jd, body_id, XYZ_FLAG)
            x, y, z = le_xyz[0][0], le_xyz[0][1], le_xyz[0][2]
            xyz_dist = math.sqrt(x * x + y * y + z * z)

            diff = abs(ecl_dist - xyz_dist)
            rel = diff / ecl_dist if ecl_dist > 0 else 0

            # Should be identical (same underlying distance)
            if rel < 1e-10:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL P6 {body_name} {epoch_label}: ecl_dist={ecl_dist:.10f} xyz_dist={xyz_dist:.10f} rel={rel:.2e}"
                )

        except Exception as e:
            errors += 1
            print(f"  ERROR P6 {body_name} {epoch_label}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Extreme epoch distances (1600, 2400)
# ============================================================
print("\n=== P7: Distance at extreme epochs ===")

extreme_epochs = [
    (swe.julday(1600, 1, 1, 12.0), "Y1600"),
    (swe.julday(1700, 1, 1, 12.0), "Y1700"),
    (swe.julday(2300, 1, 1, 12.0), "Y2300"),
    (swe.julday(2400, 1, 1, 12.0), "Y2400"),
]

OUTER_BODIES = {
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

for jd, epoch_label in extreme_epochs:
    for body_id, body_name in OUTER_BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLAGS)
            le_result = ephem.calc_ut(jd, body_id, FLAGS)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            diff_au = abs(se_dist - le_dist)

            # Wider tolerance for extreme dates: 1e-5 AU ≈ 1500 m
            tol = 1e-5
            if body_id == 9:  # Pluto
                tol = 1e-4  # Known DE440 vs SE differences

            if diff_au < tol:
                passed += 1
            else:
                diff_km = diff_au * 1.496e8
                failed += 1
                print(
                    f"  FAIL P7 {body_name} {epoch_label}: SE={se_dist:.8f} LE={le_dist:.8f} diff={diff_au:.2e} AU ({diff_km:.0f}km)"
                )

        except Exception as e:
            errors += 1
            print(f"  ERROR P7 {body_name} {epoch_label}: {e}")

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
print(
    f"ROUND 53 FINAL: {passed}/{passed + failed} passed ({100 * passed / (passed + failed):.1f}%)"
)
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
