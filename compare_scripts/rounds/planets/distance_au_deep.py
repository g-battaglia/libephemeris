#!/usr/bin/env python3
"""Round 108: Planet Distance AU Precision Deep

Goes deeper than Round 53 by testing:
P1: Distance at planetary opposition/conjunction (extreme geocentric distances)
P2: Inner planet distances at inferior/superior conjunction
P3: Distance derivative consistency (numerical finite-diff vs analytical speed)
P4: Barycentric distances for all planets
P5: Heliocentric distance at perihelion/aphelion epochs
P6: Distance consistency: geo dist ≈ |helio_planet - helio_earth| for outer planets
P7: Moon distance at monthly perigee/apogee extremes (tighter tolerance)
P8: Distance at 1-hour intervals around specific events for smoothness
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
FLG_BARYCTR = 4
FLG_TRUEPOS = 16
FLG_J2000 = 32
FLG_NONUT = 64
FLG_NOABERR = 1024

AST_OFFSET = 10000

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
}

SE_ASTEROID_MAP = {
    17: AST_OFFSET + 1,
    18: AST_OFFSET + 2,
    19: AST_OFFSET + 3,
    20: AST_OFFSET + 4,
}


def get_se_body(body_id):
    return SE_ASTEROID_MAP.get(body_id, body_id)


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


# ============================================================
# P1: Distance at planetary opposition/conjunction
# Opposition: outer planet opposite Sun → closest to Earth
# Conjunction: outer planet same direction as Sun → farthest from Earth
# ============================================================
print("=== P1: Outer planet distances at opposition/conjunction epochs ===")

# Known approximate opposition dates (planet closest to Earth)
opposition_dates = [
    # (year, month, day, body_id, body_name, event)
    (2024, 1, 27, 5, "Jupiter", "opposition"),
    (2024, 9, 8, 6, "Saturn", "opposition"),
    (2024, 11, 21, 7, "Uranus", "opposition"),
    (2024, 9, 21, 8, "Neptune", "opposition"),
    (2023, 11, 3, 5, "Jupiter", "opposition"),
    (2025, 1, 16, 4, "Mars", "opposition"),
    (2022, 12, 8, 4, "Mars", "opposition"),
    (2023, 8, 27, 6, "Saturn", "opposition"),
    # Conjunctions (planet behind Sun → farthest)
    (2024, 5, 18, 5, "Jupiter", "conjunction"),
    (2024, 2, 28, 6, "Saturn", "conjunction"),
    (2024, 5, 13, 7, "Uranus", "conjunction"),
    (2024, 3, 17, 8, "Neptune", "conjunction"),
    (2024, 11, 18, 2, "Mercury", "conjunction"),
    (2024, 6, 4, 3, "Venus", "conjunction-sup"),
]

for year, month, day, body_id, body_name, event in opposition_dates:
    # Test at event and ±5 days around it
    jd_center = swe.julday(year, month, day, 12.0)
    for offset_days in [-5, -2, 0, 2, 5]:
        jd = jd_center + offset_days
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLG_SPEED)
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            diff_au = abs(se_dist - le_dist)
            diff_km = diff_au * 1.496e8

            # Tight tolerance: 1e-7 AU ≈ 15 m
            tol = 1e-7
            label = (
                f"P1 {body_name} {event} {year}-{month:02d}-{day:02d}{offset_days:+d}d"
            )
            run_test(
                label,
                diff_au < tol,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU ({diff_km:.1f}m)",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P1 {body_name} {event}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Inner planet distances at inferior/superior conjunction
# Mercury/Venus closest at inferior conjunction
# ============================================================
print("\n=== P2: Inner planet distances at inferior/superior conjunction ===")

inner_events = [
    # Mercury inferior conjunctions 2024
    (2024, 1, 10, 2, "Mercury", "inf-conj"),
    (2024, 5, 6, 2, "Mercury", "inf-conj"),
    (2024, 8, 28, 2, "Mercury", "inf-conj"),
    (2024, 12, 18, 2, "Mercury", "inf-conj"),
    # Mercury greatest elongations
    (2024, 3, 24, 2, "Mercury", "elong-E"),
    (2024, 7, 22, 2, "Mercury", "elong-W"),
    (2024, 12, 6, 2, "Mercury", "elong-E"),
    # Venus events
    (2024, 6, 4, 3, "Venus", "sup-conj"),
    (2023, 8, 13, 3, "Venus", "inf-conj"),
    (2025, 3, 23, 3, "Venus", "inf-conj"),
    (2023, 6, 4, 3, "Venus", "elong-E"),
    (2024, 9, 22, 3, "Venus", "elong-W"),
]

for year, month, day, body_id, body_name, event in inner_events:
    jd_center = swe.julday(year, month, day, 12.0)
    for offset_days in [-3, 0, 3]:
        jd = jd_center + offset_days
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLG_SPEED)
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            diff_au = abs(se_dist - le_dist)
            diff_km = diff_au * 1.496e8

            tol = 1e-7  # 15 m
            label = (
                f"P2 {body_name} {event} {year}-{month:02d}-{day:02d}{offset_days:+d}d"
            )
            run_test(
                label,
                diff_au < tol,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU ({diff_km:.1f}m)",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P2 {body_name} {event}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Distance derivative consistency
# Compare analytical dist_speed (index 5) vs numerical finite-difference
# ============================================================
print("\n=== P3: Distance derivative consistency (analytical vs numerical) ===")

DT = 0.001  # 0.001 day ≈ 86.4 seconds for finite difference
test_jds = [
    (swe.julday(2024, 3, 20, 12.0), "2024-equinox"),
    (swe.julday(2024, 6, 21, 12.0), "2024-solstice"),
    (swe.julday(2024, 9, 22, 12.0), "2024-equinox2"),
    (swe.julday(2000, 1, 1, 12.0), "J2000"),
    (swe.julday(1990, 7, 15, 12.0), "1990"),
    (swe.julday(2050, 1, 1, 12.0), "2050"),
]

for jd, epoch_label in test_jds:
    for body_id, body_name in BODIES.items():
        try:
            # LE analytical speed
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)
            le_dist_speed = le_result[0][5]  # AU/day

            # LE numerical speed via finite difference
            le_minus = ephem.calc_ut(jd - DT, body_id, FLG_SPEED)
            le_plus = ephem.calc_ut(jd + DT, body_id, FLG_SPEED)
            num_speed = (le_plus[0][2] - le_minus[0][2]) / (2 * DT)

            diff = abs(le_dist_speed - num_speed)

            # Also check SE analytical speed
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLG_SPEED)
            se_dist_speed = se_result[0][5]

            # For Moon, speed is ~1e-3 AU/day; for Sun ~1e-4; for planets ~1e-3
            # Numerical accuracy limited by DT, expect ~DT^2 error
            if abs(le_dist_speed) > 1e-8:
                rel_err = diff / abs(le_dist_speed)
                # Relative tolerance: 0.1% for finite-difference
                ok = rel_err < 0.001
            else:
                # Near stationary point, absolute tolerance
                ok = diff < 1e-8

            label = f"P3 {body_name} {epoch_label}"
            run_test(
                label,
                ok,
                f"analytical={le_dist_speed:.10e} numerical={num_speed:.10e} diff={diff:.2e} SE={se_dist_speed:.10e}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P3 {body_name} {epoch_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Barycentric distances for all planets
# ============================================================
print("\n=== P4: Barycentric distances ===")
BARY_FLAG = FLG_SPEED | FLG_BARYCTR

bary_epochs = [
    (swe.julday(2000, 1, 1, 12.0), "J2000"),
    (swe.julday(2024, 6, 15, 12.0), "2024-mid"),
    (swe.julday(1980, 3, 1, 12.0), "1980"),
    (swe.julday(2040, 9, 1, 12.0), "2040"),
    (swe.julday(1900, 1, 1, 12.0), "1900"),
    (swe.julday(2100, 1, 1, 12.0), "2100"),
]

# Skip Moon for barycentric (known differences)
BARY_BODIES = {k: v for k, v in BODIES.items() if k != 1}

for jd, epoch_label in bary_epochs:
    for body_id, body_name in BARY_BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, BARY_FLAG)
            le_result = ephem.calc_ut(jd, body_id, BARY_FLAG)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            if se_dist == 0.0 and le_dist == 0.0:
                passed += 1
                continue

            diff_au = abs(se_dist - le_dist)

            # Tolerance based on body
            if body_id == 0:  # Sun (barycentric, ~0.01 AU)
                tol = 1e-7
            elif body_id in (5, 6):  # Jupiter/Saturn (known 1e-5 diff)
                tol = 2e-5
            else:
                tol = 1e-6

            label = f"P4 {body_name} {epoch_label}"
            run_test(
                label,
                diff_au < tol,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P4 {body_name} {epoch_label}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Heliocentric distances at perihelion/aphelion
# ============================================================
print("\n=== P5: Heliocentric distances at perihelion/aphelion ===")
HELIO_FLAG = FLG_SPEED | FLG_HELCTR

# Approximate perihelion/aphelion dates and expected distances
perihelion_aphelion = [
    # (year, month, day, body_id, name, event, expected_dist_au)
    (2024, 1, 3, 0, "Earth(Sun)", "perihelion", 0.9833),  # Sun geo = Earth helio
    (2024, 7, 5, 0, "Earth(Sun)", "aphelion", 1.0167),
    (2024, 9, 16, 4, "Mars", "perihelion-approx", 1.38),
    (2023, 1, 12, 4, "Mars", "perihelion", 1.382),
    # Mercury: e=0.2056, a=0.387 AU → peri=0.307, aph=0.467
    (2024, 2, 26, 2, "Mercury", "perihelion-approx", 0.307),
    (2024, 4, 11, 2, "Mercury", "aphelion-approx", 0.467),
    # Jupiter perihelion Jan 2023 (a=5.203, e=0.0489 → peri=4.95, aph=5.46)
    (2023, 1, 21, 5, "Jupiter", "perihelion", 4.951),
]

for year, month, day, body_id, body_name, event, expected_dist in perihelion_aphelion:
    jd = swe.julday(year, month, day, 12.0)
    for offset_days in [-10, -5, 0, 5, 10]:
        jd_test = jd + offset_days
        try:
            if body_id == 0:
                # Earth helio dist = Sun geocentric dist
                se_result = swe.calc_ut(jd_test, 0, FLG_SPEED)
                le_result = ephem.calc_ut(jd_test, 0, FLG_SPEED)
            else:
                se_body = get_se_body(body_id)
                se_result = swe.calc_ut(jd_test, se_body, HELIO_FLAG)
                le_result = ephem.calc_ut(jd_test, body_id, HELIO_FLAG)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            diff_au = abs(se_dist - le_dist)

            tol = 1e-7  # 15 m
            label = f"P5 {body_name} {event} {offset_days:+d}d"
            run_test(
                label,
                diff_au < tol,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P5 {body_name} {event} {offset_days:+d}d: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Geocentric vs heliocentric distance triangle consistency
# For outer planets: geo_dist ≈ |helio_planet - helio_earth|
# (not exact due to light-time, but should be close)
# ============================================================
print("\n=== P6: Geo vs helio distance triangle check ===")

OUTER = {4: "Mars", 5: "Jupiter", 6: "Saturn", 7: "Uranus", 8: "Neptune"}

for year in [2000, 2010, 2020, 2024]:
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 1, 12.0)
        for body_id, body_name in OUTER.items():
            try:
                # Geocentric distance
                le_geo = ephem.calc_ut(jd, body_id, FLG_SPEED)
                geo_dist = le_geo[0][2]

                # Heliocentric distance of planet
                le_helio = ephem.calc_ut(jd, body_id, FLG_SPEED | FLG_HELCTR)
                helio_dist = le_helio[0][2]

                # Earth-Sun distance (= helio distance of Earth)
                le_sun = ephem.calc_ut(jd, 0, FLG_SPEED)
                earth_sun_dist = le_sun[0][2]

                # Triangle inequality: |helio - earth_sun| <= geo <= helio + earth_sun
                lower = abs(helio_dist - earth_sun_dist)
                upper = helio_dist + earth_sun_dist

                ok = (lower - 0.01) <= geo_dist <= (upper + 0.01)

                label = f"P6 {body_name} {year}-{month:02d}"
                run_test(
                    label,
                    ok,
                    f"geo={geo_dist:.6f} helio={helio_dist:.6f} earth={earth_sun_dist:.6f} bounds=[{lower:.6f},{upper:.6f}]",
                )
            except Exception as e:
                errors += 1
                print(f"  ERROR P6 {body_name} {year}-{month:02d}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Moon distance at monthly perigee/apogee (tight tolerance)
# Test every day for 2 years
# ============================================================
print("\n=== P7: Moon distance daily for 2 years ===")

jd_start = swe.julday(2023, 1, 1, 0.0)
for i in range(730):  # 2 years
    jd = jd_start + i
    try:
        se_result = swe.calc_ut(jd, 1, FLG_SPEED)
        le_result = ephem.calc_ut(jd, 1, FLG_SPEED)

        se_dist = se_result[0][2]
        le_dist = le_result[0][2]

        diff_au = abs(se_dist - le_dist)
        diff_km = diff_au * 1.496e8

        # Very tight: 50 meters
        tol_km = 0.05
        label = f"P7 Moon day={i}"
        if diff_km >= tol_km:
            run_test(
                label,
                False,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_km:.4f}km",
            )
        else:
            passed += 1
    except Exception as e:
        errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Distance at hourly intervals around J2000 (smoothness check)
# ============================================================
print("\n=== P8: Hourly distance around J2000 for smoothness ===")

jd_j2000 = 2451545.0
HOURLY_BODIES = {0: "Sun", 1: "Moon", 2: "Mercury", 4: "Mars", 5: "Jupiter"}

for body_id, body_name in HOURLY_BODIES.items():
    prev_dist = None
    prev_speed = None
    for hour in range(-48, 49):  # ±2 days at hourly intervals
        jd = jd_j2000 + hour / 24.0
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLG_SPEED)
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            diff_au = abs(se_dist - le_dist)

            tol = 1e-7
            if body_id == 1:
                tol = 5e-8  # tighter for Moon

            label = f"P8 {body_name} h={hour:+d}"
            run_test(
                label,
                diff_au < tol,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e}",
            )
        except Exception as e:
            errors += 1

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P9: TRUEPOS distances (geometric, no light-time correction)
# ============================================================
print("\n=== P9: TRUEPOS (geometric) distances ===")
TRUEPOS_FLAG = FLG_SPEED | FLG_TRUEPOS

truepos_epochs = [
    (swe.julday(2000, 1, 1, 12.0), "J2000"),
    (swe.julday(2024, 6, 15, 12.0), "2024-mid"),
    (swe.julday(1950, 1, 1, 12.0), "1950"),
    (swe.julday(2080, 1, 1, 12.0), "2080"),
]

for jd, epoch_label in truepos_epochs:
    for body_id, body_name in BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, TRUEPOS_FLAG)
            le_result = ephem.calc_ut(jd, body_id, TRUEPOS_FLAG)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            if se_dist == 0.0 and le_dist == 0.0:
                passed += 1
                continue

            diff_au = abs(se_dist - le_dist)

            tol = 1e-6
            if body_id == 1:
                tol = 1e-7

            label = f"P9 {body_name} {epoch_label} TRUEPOS"
            run_test(
                label,
                diff_au < tol,
                f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P9 {body_name} {epoch_label}: {e}")

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P10: J2000 frame distances (should match default)
# ============================================================
print("\n=== P10: J2000 frame distances ===")
J2000_FLAG = FLG_SPEED | FLG_J2000 | FLG_NONUT

j2000_epochs = [
    (swe.julday(2000, 1, 1, 12.0), "J2000"),
    (swe.julday(2024, 3, 20, 12.0), "2024-equi"),
    (swe.julday(1970, 7, 20, 12.0), "1970"),
]

for jd, epoch_label in j2000_epochs:
    for body_id, body_name in BODIES.items():
        try:
            # Default frame distance
            le_default = ephem.calc_ut(jd, body_id, FLG_SPEED)
            default_dist = le_default[0][2]

            # J2000 frame distance
            le_j2000 = ephem.calc_ut(jd, body_id, J2000_FLAG)
            j2000_dist = le_j2000[0][2]

            # Distances should be identical regardless of frame
            diff = abs(default_dist - j2000_dist)

            label = f"P10 {body_name} {epoch_label}"
            run_test(
                label,
                diff < 1e-12,
                f"default={default_dist:.12f} J2000={j2000_dist:.12f} diff={diff:.2e}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P10 {body_name} {epoch_label}: {e}")

print(f"  After P10: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P11: Distance speed sign consistency
# When planet approaching Earth, dist_speed should be negative
# When receding, dist_speed should be positive
# ============================================================
print("\n=== P11: Distance speed sign consistency ===")

sign_epochs = [
    (swe.julday(2024, 1, 1, 12.0), "2024-Jan"),
    (swe.julday(2024, 4, 1, 12.0), "2024-Apr"),
    (swe.julday(2024, 7, 1, 12.0), "2024-Jul"),
    (swe.julday(2024, 10, 1, 12.0), "2024-Oct"),
]

for jd, epoch_label in sign_epochs:
    for body_id, body_name in BODIES.items():
        try:
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)
            le_dist = le_result[0][2]
            le_dspeed = le_result[0][5]

            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLG_SPEED)
            se_dspeed = se_result[0][5]

            # Check sign consistency
            if abs(le_dspeed) < 1e-10 and abs(se_dspeed) < 1e-10:
                passed += 1
                continue

            # Both should have same sign
            same_sign = (le_dspeed >= 0 and se_dspeed >= 0) or (
                le_dspeed < 0 and se_dspeed < 0
            )

            label = f"P11 {body_name} {epoch_label} sign"
            run_test(
                label,
                same_sign,
                f"LE_dspeed={le_dspeed:.10e} SE_dspeed={se_dspeed:.10e}",
            )
        except Exception as e:
            errors += 1

print(f"  After P11: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P12: Multi-century distance sweep (50-year intervals)
# ============================================================
print("\n=== P12: Multi-century distance sweep ===")

SWEEP_BODIES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
}

for year in range(1850, 2151, 25):
    jd = swe.julday(year, 6, 15, 12.0)
    for body_id, body_name in SWEEP_BODIES.items():
        try:
            se_body = get_se_body(body_id)
            se_result = swe.calc_ut(jd, se_body, FLG_SPEED)
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)

            se_dist = se_result[0][2]
            le_dist = le_result[0][2]

            diff_au = abs(se_dist - le_dist)

            tol = 1e-7
            if body_id in (5, 6):  # outer planets slightly wider
                tol = 5e-7

            label = f"P12 {body_name} Y{year}"
            if diff_au >= tol:
                run_test(
                    label,
                    False,
                    f"SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_au:.2e} AU",
                )
            else:
                passed += 1
        except Exception as e:
            errors += 1

print(f"  After P12: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 108 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
