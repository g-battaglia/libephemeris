#!/usr/bin/env python3
"""Round 111: Planetary Magnitude/Phase Deep

Deep test of pheno_ut for all planets across multiple epochs.
P1: Phase angle for all planets at multiple epochs
P2: Elongation for all planets
P3: Apparent diameter for all planets
P4: Magnitude for all planets
P5: Phase (illuminated fraction) for all planets
P6: Moon phase across full lunation cycle
P7: Mercury/Venus phase near inferior conjunction
P8: Heliocentric pheno (elongation should be 0)
P9: Phase consistency: phase_angle and elongation relationships
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

FLG_SPEED = 256
FLG_HELCTR = 8
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


# Test epochs
TEST_EPOCHS = [
    (swe.julday(2000, 1, 1, 12.0), "J2000"),
    (swe.julday(2024, 1, 15, 12.0), "2024-Jan"),
    (swe.julday(2024, 4, 15, 12.0), "2024-Apr"),
    (swe.julday(2024, 7, 15, 12.0), "2024-Jul"),
    (swe.julday(2024, 10, 15, 12.0), "2024-Oct"),
    (swe.julday(1990, 6, 15, 12.0), "1990"),
    (swe.julday(2050, 3, 15, 12.0), "2050"),
    (swe.julday(1950, 9, 15, 12.0), "1950"),
]


def get_pheno(jd, body_id):
    """Get pheno values from both SE and LE."""
    se_body = get_se_body(body_id)

    # SE: returns flat tuple of 20 values
    se_raw = swe.pheno_ut(jd, se_body, FLG_SPEED)
    se_phase_angle = se_raw[0]
    se_phase = se_raw[1]
    se_elongation = se_raw[2]
    se_diameter = se_raw[3]
    se_magnitude = se_raw[4]

    # LE: returns ((phase_angle, phase, elongation, diameter, magnitude), retflag)
    le_result = ephem.pheno_ut(jd, body_id, FLG_SPEED)
    le_vals = le_result[0]
    le_phase_angle = le_vals[0]
    le_phase = le_vals[1]
    le_elongation = le_vals[2]
    le_diameter = le_vals[3]
    le_magnitude = le_vals[4]

    return (
        (se_phase_angle, se_phase, se_elongation, se_diameter, se_magnitude),
        (le_phase_angle, le_phase, le_elongation, le_diameter, le_magnitude),
    )


# ============================================================
# P1: Phase angle for all planets
# ============================================================
print("=== P1: Phase angle for all planets ===")

for jd, epoch_label in TEST_EPOCHS:
    for body_id, body_name in BODIES.items():
        try:
            se_vals, le_vals = get_pheno(jd, body_id)
            se_pa, le_pa = se_vals[0], le_vals[0]

            diff = abs(se_pa - le_pa)
            # Phase angle in degrees, tolerance 0.01° (36")
            tol = 0.01
            if body_id in (7, 8, 9):  # outer planets — wider
                tol = 0.05

            label = f"P1 {body_name} {epoch_label} phase_angle"
            run_test(
                label,
                diff < tol,
                f"SE={se_pa:.6f}° LE={le_pa:.6f}° diff={diff:.6f}°",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P1 {body_name} {epoch_label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Elongation for all planets
# ============================================================
print("\n=== P2: Elongation for all planets ===")

for jd, epoch_label in TEST_EPOCHS:
    for body_id, body_name in BODIES.items():
        if body_id == 0:  # Sun elongation is 0 by definition
            continue
        try:
            se_vals, le_vals = get_pheno(jd, body_id)
            se_el, le_el = se_vals[2], le_vals[2]

            diff = abs(se_el - le_el)
            tol = 0.01  # 36"

            label = f"P2 {body_name} {epoch_label} elongation"
            run_test(
                label,
                diff < tol,
                f"SE={se_el:.6f}° LE={le_el:.6f}° diff={diff:.6f}°",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P2 {body_name} {epoch_label}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Apparent diameter for all planets
# ============================================================
print("\n=== P3: Apparent diameter for all planets ===")

for jd, epoch_label in TEST_EPOCHS:
    for body_id, body_name in BODIES.items():
        try:
            se_vals, le_vals = get_pheno(jd, body_id)
            se_diam, le_diam = se_vals[3], le_vals[3]

            if se_diam == 0.0 and le_diam == 0.0:
                passed += 1
                continue

            if se_diam == 0.0 or le_diam == 0.0:
                run_test(
                    f"P3 {body_name} {epoch_label} diam",
                    False,
                    f"SE={se_diam:.6f} LE={le_diam:.6f} (one is zero)",
                )
                continue

            rel_diff = abs(se_diam - le_diam) / se_diam

            # Diameter tolerance: 1% relative
            tol = 0.01

            label = f"P3 {body_name} {epoch_label} diameter"
            run_test(
                label,
                rel_diff < tol,
                f"SE={se_diam:.6f}° LE={le_diam:.6f}° rel={rel_diff:.4f}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P3 {body_name} {epoch_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Magnitude for all planets
# ============================================================
print("\n=== P4: Magnitude for all planets ===")

for jd, epoch_label in TEST_EPOCHS:
    for body_id, body_name in BODIES.items():
        if body_id == 0:  # Sun magnitude is special
            continue
        try:
            se_vals, le_vals = get_pheno(jd, body_id)
            se_mag, le_mag = se_vals[4], le_vals[4]

            diff = abs(se_mag - le_mag)

            # Magnitude tolerance: 0.2 mag (generous for model differences)
            tol = 0.2
            if body_id == 1:  # Moon — model-dependent
                tol = 0.5
            if body_id in (7, 8, 9):  # outer — model-dependent
                tol = 0.5

            label = f"P4 {body_name} {epoch_label} magnitude"
            run_test(
                label,
                diff < tol,
                f"SE={se_mag:.3f} LE={le_mag:.3f} diff={diff:.3f}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P4 {body_name} {epoch_label}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Phase (illuminated fraction)
# ============================================================
print("\n=== P5: Phase (illuminated fraction) ===")

for jd, epoch_label in TEST_EPOCHS:
    for body_id, body_name in BODIES.items():
        try:
            se_vals, le_vals = get_pheno(jd, body_id)
            se_phase, le_phase = se_vals[1], le_vals[1]

            diff = abs(se_phase - le_phase)

            # Phase tolerance: 0.01 (1%)
            tol = 0.01
            if body_id == 0:  # Sun phase is 0.0
                tol = 0.001

            label = f"P5 {body_name} {epoch_label} phase"
            run_test(
                label,
                diff < tol,
                f"SE={se_phase:.6f} LE={le_phase:.6f} diff={diff:.6f}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P5 {body_name} {epoch_label}: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Moon phase across full lunation cycle (29.53 days)
# ============================================================
print("\n=== P6: Moon phase across lunation cycle ===")

# Start near New Moon January 2024
jd_new_moon = swe.julday(2024, 1, 11, 12.0)

for i in range(60):  # 2 lunations at half-day intervals
    jd = jd_new_moon + i * 0.5
    try:
        se_vals, le_vals = get_pheno(jd, 1)

        # Phase angle
        diff_pa = abs(se_vals[0] - le_vals[0])
        run_test(
            f"P6 Moon lunation d={i * 0.5:.1f} PA",
            diff_pa < 0.01,
            f"SE={se_vals[0]:.4f} LE={le_vals[0]:.4f} diff={diff_pa:.4f}",
        )

        # Phase (illuminated fraction)
        diff_ph = abs(se_vals[1] - le_vals[1])
        run_test(
            f"P6 Moon lunation d={i * 0.5:.1f} phase",
            diff_ph < 0.01,
            f"SE={se_vals[1]:.6f} LE={le_vals[1]:.6f} diff={diff_ph:.6f}",
        )

        # Elongation
        diff_el = abs(se_vals[2] - le_vals[2])
        run_test(
            f"P6 Moon lunation d={i * 0.5:.1f} elong",
            diff_el < 0.01,
            f"SE={se_vals[2]:.4f} LE={le_vals[2]:.4f} diff={diff_el:.4f}",
        )
    except Exception as e:
        errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Mercury/Venus phase near inferior conjunction
# ============================================================
print("\n=== P7: Mercury/Venus near inferior conjunction ===")

# Mercury inferior conjunctions in 2024
mercury_inf_conj = [
    swe.julday(2024, 1, 10, 0.0),
    swe.julday(2024, 5, 6, 0.0),
    swe.julday(2024, 8, 28, 0.0),
    swe.julday(2024, 12, 18, 0.0),
]

for jd_conj in mercury_inf_conj:
    for offset in range(-10, 11, 2):
        jd = jd_conj + offset
        try:
            se_vals, le_vals = get_pheno(jd, 2)

            # Phase angle
            diff_pa = abs(se_vals[0] - le_vals[0])
            tol_pa = 0.05  # wider near conjunction
            if diff_pa >= tol_pa:
                run_test(
                    f"P7 Mercury inf_conj d={offset:+d} PA",
                    False,
                    f"SE={se_vals[0]:.4f} LE={le_vals[0]:.4f}",
                )
            else:
                passed += 1

            # Elongation
            diff_el = abs(se_vals[2] - le_vals[2])
            if diff_el >= 0.05:
                run_test(
                    f"P7 Mercury inf_conj d={offset:+d} elong",
                    False,
                    f"SE={se_vals[2]:.4f} LE={le_vals[2]:.4f}",
                )
            else:
                passed += 1
        except Exception as e:
            errors += 1

# Venus inferior conjunction (Aug 2023)
jd_venus_inf = swe.julday(2023, 8, 13, 0.0)
for offset in range(-15, 16, 3):
    jd = jd_venus_inf + offset
    try:
        se_vals, le_vals = get_pheno(jd, 3)

        diff_pa = abs(se_vals[0] - le_vals[0])
        if diff_pa >= 0.05:
            run_test(
                f"P7 Venus inf_conj d={offset:+d} PA",
                False,
                f"SE={se_vals[0]:.4f} LE={le_vals[0]:.4f}",
            )
        else:
            passed += 1

        diff_el = abs(se_vals[2] - le_vals[2])
        if diff_el >= 0.05:
            run_test(
                f"P7 Venus inf_conj d={offset:+d} elong",
                False,
                f"SE={se_vals[2]:.4f} LE={le_vals[2]:.4f}",
            )
        else:
            passed += 1
    except Exception as e:
        errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Phase sanity checks
# ============================================================
print("\n=== P8: Phase sanity checks ===")

for jd, epoch_label in TEST_EPOCHS[:4]:
    for body_id, body_name in BODIES.items():
        try:
            le_result = ephem.pheno_ut(jd, body_id, FLG_SPEED)
            le_vals = le_result[0]

            pa = le_vals[0]  # phase angle
            phase = le_vals[1]  # illuminated fraction
            elong = le_vals[2]  # elongation
            diam = le_vals[3]  # diameter
            mag = le_vals[4]  # magnitude

            # Phase angle should be 0-180
            run_test(
                f"P8 {body_name} {epoch_label} PA_range",
                0 <= pa <= 180,
                f"PA={pa:.4f} out of [0,180]",
            )

            # Phase should be 0-1 (or close)
            if body_id != 0:  # Sun phase is special
                run_test(
                    f"P8 {body_name} {epoch_label} phase_range",
                    -0.01 <= phase <= 1.01,
                    f"phase={phase:.6f} out of [0,1]",
                )

            # Elongation should be 0-180
            if body_id != 0:
                run_test(
                    f"P8 {body_name} {epoch_label} elong_range",
                    0 <= elong <= 180.01,
                    f"elong={elong:.4f} out of [0,180]",
                )

            # Diameter should be non-negative
            run_test(
                f"P8 {body_name} {epoch_label} diam_nonneg",
                diam >= 0,
                f"diam={diam:.6f} negative",
            )

        except Exception as e:
            errors += 1
            print(f"  ERROR P8 {body_name} {epoch_label}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P9: Multi-year sweep for outer planet phenomena
# ============================================================
print("\n=== P9: Multi-year outer planet phenomena ===")

OUTER = {4: "Mars", 5: "Jupiter", 6: "Saturn"}

for year in range(2000, 2026):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 1, 12.0)
        for body_id, body_name in OUTER.items():
            try:
                se_vals, le_vals = get_pheno(jd, body_id)

                # Phase angle
                diff_pa = abs(se_vals[0] - le_vals[0])
                if diff_pa >= 0.02:
                    run_test(
                        f"P9 {body_name} {year}-{month:02d} PA",
                        False,
                        f"SE={se_vals[0]:.4f} LE={le_vals[0]:.4f} diff={diff_pa:.4f}",
                    )
                else:
                    passed += 1

                # Magnitude
                diff_mag = abs(se_vals[4] - le_vals[4])
                tol_mag = 0.15
                if diff_mag >= tol_mag:
                    run_test(
                        f"P9 {body_name} {year}-{month:02d} mag",
                        False,
                        f"SE={se_vals[4]:.3f} LE={le_vals[4]:.3f} diff={diff_mag:.3f}",
                    )
                else:
                    passed += 1

            except Exception as e:
                errors += 1

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P10: Sun diameter consistency across year
# ============================================================
print("\n=== P10: Sun diameter across year ===")

jd_start = swe.julday(2024, 1, 1, 12.0)
for day in range(0, 365, 5):
    jd = jd_start + day
    try:
        se_vals, le_vals = get_pheno(jd, 0)
        se_diam, le_diam = se_vals[3], le_vals[3]

        rel_diff = abs(se_diam - le_diam) / se_diam if se_diam > 0 else 0

        if rel_diff >= 0.005:  # 0.5%
            run_test(
                f"P10 Sun diam day={day}",
                False,
                f"SE={se_diam:.8f} LE={le_diam:.8f} rel={rel_diff:.6f}",
            )
        else:
            passed += 1
    except Exception as e:
        errors += 1

print(f"  After P10: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 111 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
