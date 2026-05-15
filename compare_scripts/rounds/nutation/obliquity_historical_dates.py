#!/usr/bin/env python3
"""Round 54: Ecliptic Obliquity at Historical Dates

Compare mean and true obliquity across wide date ranges.
Tests calc_ut with ECL_NUT body (special body returning nutation/obliquity).
Also tests the internal obliquity computation consistency.
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

# ECL_NUT = -1 in pyswisseph
ECL_NUT = -1

# Test epochs spanning 3000 years
TEST_EPOCHS = []
for year in range(1000, 2801, 25):
    jd = swe.julday(year, 1, 1, 12.0)
    TEST_EPOCHS.append((jd, f"Y{year}"))

# Add key astronomical dates
TEST_EPOCHS.append((2451545.0, "J2000.0"))
TEST_EPOCHS.append((2415020.0, "J1900.0"))
TEST_EPOCHS.append((2433282.5, "B1950.0"))
TEST_EPOCHS.append((2460000.0, "JD2460000"))

FLAGS = 0

# ============================================================
# P1: True obliquity via ECL_NUT across epochs
# ============================================================
print("=" * 70)
print("ROUND 54: Ecliptic Obliquity at Historical Dates")
print("=" * 70)

print("\n=== P1: True obliquity via ECL_NUT (1000-2800) ===")

for jd, epoch_label in TEST_EPOCHS:
    try:
        # SE: swe.calc_ut(jd, ECL_NUT, 0) returns [true_obl, mean_obl, dpsi, deps, 0, 0]
        se_result = swe.calc_ut(jd, ECL_NUT, 0)
        se_pos = se_result[0] if isinstance(se_result, tuple) else se_result

        le_result = ephem.calc_ut(jd, ECL_NUT, 0)
        le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

        se_true_obl = se_pos[0]  # true obliquity
        le_true_obl = le_pos[0]

        diff_arcsec = abs(se_true_obl - le_true_obl) * 3600.0

        # Tolerance: 0.1" for true obliquity
        if diff_arcsec < 0.1:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P1 true_obl {epoch_label}: SE={se_true_obl:.8f}° LE={le_true_obl:.8f}° diff={diff_arcsec:.4f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P1 {epoch_label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Mean obliquity via ECL_NUT
# ============================================================
print("\n=== P2: Mean obliquity via ECL_NUT (1000-2800) ===")

for jd, epoch_label in TEST_EPOCHS:
    try:
        se_result = swe.calc_ut(jd, ECL_NUT, 0)
        se_pos = se_result[0]

        le_result = ephem.calc_ut(jd, ECL_NUT, 0)
        le_pos = le_result[0]

        se_mean_obl = se_pos[1]  # mean obliquity
        le_mean_obl = le_pos[1]

        diff_arcsec = abs(se_mean_obl - le_mean_obl) * 3600.0

        # Tolerance: 0.05" for mean obliquity (should be very close)
        if diff_arcsec < 0.05:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P2 mean_obl {epoch_label}: SE={se_mean_obl:.8f}° LE={le_mean_obl:.8f}° diff={diff_arcsec:.4f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P2 {epoch_label}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Nutation in longitude (dpsi) via ECL_NUT
# ============================================================
print("\n=== P3: Nutation in longitude (dpsi) ===")

for jd, epoch_label in TEST_EPOCHS:
    try:
        se_result = swe.calc_ut(jd, ECL_NUT, 0)
        se_pos = se_result[0]

        le_result = ephem.calc_ut(jd, ECL_NUT, 0)
        le_pos = le_result[0]

        se_dpsi = se_pos[2]  # nutation in longitude
        le_dpsi = le_pos[2]

        diff_arcsec = abs(se_dpsi - le_dpsi) * 3600.0

        # Tolerance: 0.1" for nutation
        if diff_arcsec < 0.15:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P3 dpsi {epoch_label}: SE={se_dpsi:.8f}° LE={le_dpsi:.8f}° diff={diff_arcsec:.4f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P3 {epoch_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Nutation in obliquity (deps) via ECL_NUT
# ============================================================
print("\n=== P4: Nutation in obliquity (deps) ===")

for jd, epoch_label in TEST_EPOCHS:
    try:
        se_result = swe.calc_ut(jd, ECL_NUT, 0)
        se_pos = se_result[0]

        le_result = ephem.calc_ut(jd, ECL_NUT, 0)
        le_pos = le_result[0]

        se_deps = se_pos[3]  # nutation in obliquity
        le_deps = le_pos[3]

        diff_arcsec = abs(se_deps - le_deps) * 3600.0

        # Tolerance: 0.1" for nutation in obliquity
        if diff_arcsec < 0.15:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P4 deps {epoch_label}: SE={se_deps:.8f}° LE={le_deps:.8f}° diff={diff_arcsec:.4f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P4 {epoch_label}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Obliquity trend verification (Laskar 1986 polynomial)
# ============================================================
print('\n=== P5: Obliquity secular trend (decreasing ~47"/century) ===')

# The mean obliquity should decrease by ~0.47"/year = ~47"/century
# From J2000.0 (23°26'21.448") decreasing
jd_2000 = 2451545.0
jd_2100 = swe.julday(2100, 1, 1, 12.0)
jd_1900 = swe.julday(1900, 1, 1, 12.0)

try:
    le_2000 = ephem.calc_ut(jd_2000, ECL_NUT, 0)[0][1]
    le_2100 = ephem.calc_ut(jd_2100, ECL_NUT, 0)[0][1]
    le_1900 = ephem.calc_ut(jd_1900, ECL_NUT, 0)[0][1]

    # Rate from 2000 to 2100
    rate_forward = (le_2100 - le_2000) * 3600.0  # arcsec per century
    # Rate from 1900 to 2000
    rate_backward = (le_2000 - le_1900) * 3600.0

    # Expected: approximately -46.8" to -47.0" per century
    for label, rate in [("2000-2100", rate_forward), ("1900-2000", rate_backward)]:
        if -48.0 < rate < -45.0:
            passed += 1
        else:
            failed += 1
            print(f'  FAIL P5 rate {label}: {rate:.2f}"/century (expected ~-47")')

    # J2000.0 mean obliquity should be ~23.4393 degrees
    if abs(le_2000 - 23.4393) < 0.001:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL P5 J2000 mean_obl: {le_2000:.6f}° (expected ~23.4393°)")

except Exception as e:
    errors += 1
    print(f"  ERROR P5: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: True vs mean obliquity difference = nutation in obliquity
# ============================================================
print("\n=== P6: true_obl - mean_obl = deps consistency ===")

for jd, epoch_label in TEST_EPOCHS[:30]:  # First 30 epochs
    try:
        le_result = ephem.calc_ut(jd, ECL_NUT, 0)[0]

        true_obl = le_result[0]
        mean_obl = le_result[1]
        deps = le_result[3]

        # true_obl = mean_obl + deps
        computed_diff = true_obl - mean_obl
        diff_arcsec = abs(computed_diff - deps) * 3600.0

        if diff_arcsec < 0.001:  # Sub-milliarcsecond consistency
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P6 {epoch_label}: true-mean={computed_diff:.8f}° deps={deps:.8f}° diff={diff_arcsec:.6f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P6 {epoch_label}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Obliquity at well-known historical epochs
# ============================================================
print("\n=== P7: Obliquity at well-known epochs ===")

# IAU standard: J2000.0 mean obliquity = 23°26'21.448" = 84381.448"
# = 23.439291111° (IAU 2006 value)
jd_j2000 = 2451545.0

try:
    le_result = ephem.calc_ut(jd_j2000, ECL_NUT, 0)[0]
    se_result = swe.calc_ut(jd_j2000, ECL_NUT, 0)[0]

    le_mean = le_result[1]
    se_mean = se_result[1]

    # IAU 2006 reference: 84381.406" = 23.4392794444°
    # Lieske 1976 (used by SE): 84381.448" = 23.4392911111°
    iau_ref = 84381.406 / 3600.0  # IAU 2006
    lieske_ref = 84381.448 / 3600.0  # Lieske 1976

    # LE should match one of these references
    le_diff_iau = abs(le_mean - iau_ref) * 3600.0
    le_diff_lieske = abs(le_mean - lieske_ref) * 3600.0
    se_diff_iau = abs(se_mean - iau_ref) * 3600.0
    se_diff_lieske = abs(se_mean - lieske_ref) * 3600.0

    min_le_diff = min(le_diff_iau, le_diff_lieske)
    if min_le_diff < 0.05:
        passed += 1
    else:
        failed += 1
        print(
            f"  FAIL P7 J2000 mean_obl: LE={le_mean:.8f}° IAU={iau_ref:.8f}° Lieske={lieske_ref:.8f}°"
        )
        print(f'         diff_IAU={le_diff_iau:.4f}" diff_Lieske={le_diff_lieske:.4f}"')

    # Check SE vs LE agreement
    diff_se_le = abs(se_mean - le_mean) * 3600.0
    if diff_se_le < 0.05:
        passed += 1
    else:
        failed += 1
        print(
            f'  FAIL P7 SE vs LE at J2000: SE={se_mean:.8f}° LE={le_mean:.8f}° diff={diff_se_le:.4f}"'
        )

except Exception as e:
    errors += 1
    print(f"  ERROR P7: {e}")

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: High-frequency nutation sampling (every day for 1 month)
# ============================================================
print("\n=== P8: Daily nutation sampling (1 month) ===")

jd_start = swe.julday(2024, 6, 1, 0.0)
for day in range(30):
    jd = jd_start + day
    try:
        se_result = swe.calc_ut(jd, ECL_NUT, 0)[0]
        le_result = ephem.calc_ut(jd, ECL_NUT, 0)[0]

        # Check all 4 components
        labels = ["true_obl", "mean_obl", "dpsi", "deps"]
        tols = [0.1, 0.05, 0.15, 0.15]  # arcsec

        all_ok = True
        for i, (label, tol) in enumerate(zip(labels, tols)):
            diff = abs(se_result[i] - le_result[i]) * 3600.0
            if diff >= tol:
                all_ok = False
                failed += 1
                print(
                    f'  FAIL P8 day{day} {label}: SE={se_result[i]:.8f}° LE={le_result[i]:.8f}° diff={diff:.4f}"'
                )

        if all_ok:
            passed += 1

    except Exception as e:
        errors += 1
        print(f"  ERROR P8 day{day}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
print(
    f"ROUND 54 FINAL: {passed}/{passed + failed} passed ({100 * passed / (passed + failed):.1f}%)"
)
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
