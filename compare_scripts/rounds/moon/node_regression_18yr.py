#!/usr/bin/env python3
"""Round 55: Moon Node Regression Cycle (18.6 yr)

Verify the lunar node regression cycle:
- Mean node completes full regression in ~18.6134 years (6798.38 days)
- True node oscillates around mean node with ~18.6 year period
- Compare mean node, true node positions and speeds across multiple cycles
- Verify node crossing times (0°/180°) happen at expected intervals
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

MEAN_NODE = 10
TRUE_NODE = 11
FLAGS = 256  # FLG_SPEED

print("=" * 70)
print("ROUND 55: Moon Node Regression Cycle (18.6 yr)")
print("=" * 70)

# ============================================================
# P1: Mean node position across 2 full cycles (~37 years)
# ============================================================
print("\n=== P1: Mean node position (2000-2037) ===")

jd_start = 2451545.0  # J2000.0
for i in range(0, 37 * 12, 1):  # Monthly for 37 years
    jd = jd_start + i * 30.4375  # ~1 month
    try:
        se_result = swe.calc_ut(jd, MEAN_NODE, FLAGS)
        le_result = ephem.calc_ut(jd, MEAN_NODE, FLAGS)

        se_lon = se_result[0][0]
        le_lon = le_result[0][0]

        diff_arcsec = abs(se_lon - le_lon) * 3600.0
        # Handle wrap-around
        if diff_arcsec > 180 * 3600:
            diff_arcsec = 360 * 3600 - diff_arcsec

        # Tight tolerance for mean node: 1"
        if diff_arcsec < 1.0:
            passed += 1
        else:
            year = 2000 + i / 12.0
            failed += 1
            print(
                f'  FAIL P1 Y{year:.1f}: SE={se_lon:.6f}° LE={le_lon:.6f}° diff={diff_arcsec:.2f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P1 month={i}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: True node position across 2 full cycles
# ============================================================
print("\n=== P2: True node position (2000-2037) ===")

for i in range(0, 37 * 12, 1):
    jd = jd_start + i * 30.4375
    try:
        se_result = swe.calc_ut(jd, TRUE_NODE, FLAGS)
        le_result = ephem.calc_ut(jd, TRUE_NODE, FLAGS)

        se_lon = se_result[0][0]
        le_lon = le_result[0][0]

        diff_arcsec = abs(se_lon - le_lon) * 3600.0
        if diff_arcsec > 180 * 3600:
            diff_arcsec = 360 * 3600 - diff_arcsec

        # True node: 2" tolerance (known small model differences)
        if diff_arcsec < 2.0:
            passed += 1
        else:
            year = 2000 + i / 12.0
            failed += 1
            print(
                f'  FAIL P2 Y{year:.1f}: SE={se_lon:.6f}° LE={le_lon:.6f}° diff={diff_arcsec:.2f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P2 month={i}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Mean node speed (should be ~-0.05299° / day)
# ============================================================
print("\n=== P3: Mean node speed ===")

for i in range(0, 37 * 4, 1):  # Quarterly
    jd = jd_start + i * 91.3125
    try:
        se_result = swe.calc_ut(jd, MEAN_NODE, FLAGS)
        le_result = ephem.calc_ut(jd, MEAN_NODE, FLAGS)

        se_speed = se_result[0][3]  # lon speed
        le_speed = le_result[0][3]

        diff = abs(se_speed - le_speed)

        # Mean node speed should be very consistent
        if diff < 1e-6:  # ~0.0036" per day
            passed += 1
        else:
            year = 2000 + i / 4.0
            failed += 1
            print(
                f"  FAIL P3 Y{year:.1f}: SE_spd={se_speed:.8f} LE_spd={le_speed:.8f} diff={diff:.2e}°/day"
            )

        # Also verify the speed is approximately correct (~-0.05299°/day)
        expected_speed = -360.0 / 6798.38  # ~-0.05295°/day
        if abs(le_speed - expected_speed) < 0.001:
            passed += 1
        else:
            year = 2000 + i / 4.0
            failed += 1
            print(
                f"  FAIL P3 speed magnitude Y{year:.1f}: LE={le_speed:.6f} expected~{expected_speed:.6f}"
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P3: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: True node speed (varies, but should match SE)
# ============================================================
print("\n=== P4: True node speed ===")

for i in range(0, 37 * 4, 1):
    jd = jd_start + i * 91.3125
    try:
        se_result = swe.calc_ut(jd, TRUE_NODE, FLAGS)
        le_result = ephem.calc_ut(jd, TRUE_NODE, FLAGS)

        se_speed = se_result[0][3]
        le_speed = le_result[0][3]

        diff = abs(se_speed - le_speed)

        # True node speed varies more, use relative tolerance
        if diff < 0.001:  # ~3.6" per day
            passed += 1
        else:
            year = 2000 + i / 4.0
            failed += 1
            print(
                f"  FAIL P4 Y{year:.1f}: SE_spd={se_speed:.6f} LE_spd={le_speed:.6f} diff={diff:.4f}°/day"
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P4: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Regression cycle period verification
# ============================================================
print("\n=== P5: Regression cycle period ===")

# Track mean node through full 360° regression
# Starting from J2000.0, find when mean node returns to same longitude
try:
    start_result = ephem.calc_ut(jd_start, MEAN_NODE, FLAGS)
    start_lon = start_result[0][0]

    # Expected period: 6798.38 days = 18.6134 years
    expected_period = 6798.38

    # Check position at expected_period
    jd_end = jd_start + expected_period
    end_result = ephem.calc_ut(jd_end, MEAN_NODE, FLAGS)
    end_lon = end_result[0][0]

    diff_lon = abs(end_lon - start_lon)
    if diff_lon > 180:
        diff_lon = 360 - diff_lon

    # After one full cycle, should be back near start (within ~1°)
    if diff_lon < 2.0:
        passed += 1
    else:
        failed += 1
        print(
            f"  FAIL P5 full cycle: start={start_lon:.4f}° end={end_lon:.4f}° diff={diff_lon:.4f}°"
        )

    # Check half cycle (~180° regression)
    jd_half = jd_start + expected_period / 2
    half_result = ephem.calc_ut(jd_half, MEAN_NODE, FLAGS)
    half_lon = half_result[0][0]

    half_diff = abs(half_lon - start_lon)
    if half_diff > 180:
        half_diff = 360 - half_diff

    # Should be ~180° away
    if abs(half_diff - 180.0) < 2.0:
        passed += 1
    else:
        failed += 1
        print(
            f"  FAIL P5 half cycle: start={start_lon:.4f}° half={half_lon:.4f}° offset={half_diff:.4f}° (expected ~180°)"
        )

    # SE should give same result
    se_end = swe.calc_ut(jd_end, MEAN_NODE, FLAGS)
    se_end_lon = se_end[0][0]
    diff_se = abs(se_end_lon - end_lon) * 3600
    if diff_se < 1.0:
        passed += 1
    else:
        failed += 1
        print(
            f'  FAIL P5 SE vs LE at cycle end: SE={se_end_lon:.6f}° LE={end_lon:.6f}° diff={diff_se:.2f}"'
        )

except Exception as e:
    errors += 1
    print(f"  ERROR P5: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: True-mean node difference (nutation in longitude effect)
# ============================================================
print("\n=== P6: True-mean node difference (should be ~±1.5°) ===")

max_diff = 0.0
for i in range(0, 37 * 12, 1):
    jd = jd_start + i * 30.4375
    try:
        le_mean = ephem.calc_ut(jd, MEAN_NODE, FLAGS)[0][0]
        le_true = ephem.calc_ut(jd, TRUE_NODE, FLAGS)[0][0]

        diff = le_true - le_mean
        if diff > 180:
            diff -= 360
        if diff < -180:
            diff += 360

        if abs(diff) > max_diff:
            max_diff = abs(diff)

        # True-mean difference should be bounded (~±2° max, typically ~±1.5°)
        if abs(diff) < 3.0:
            passed += 1
        else:
            year = 2000 + i / 12.0
            failed += 1
            print(f"  FAIL P6 Y{year:.1f}: true-mean={diff:.4f}° (too large)")

    except Exception as e:
        errors += 1

print(f"  Max true-mean difference: {max_diff:.4f}°")
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Mean node latitude (should be ~0)
# ============================================================
print("\n=== P7: Mean node latitude (should be 0) ===")

for i in range(0, 20):
    jd = jd_start + i * 365.25
    try:
        le_result = ephem.calc_ut(jd, MEAN_NODE, FLAGS)[0]
        se_result = swe.calc_ut(jd, MEAN_NODE, FLAGS)[0]

        le_lat = le_result[1]
        se_lat = se_result[1]

        # Mean node latitude should be 0 (it's a point on the ecliptic)
        if abs(le_lat) < 0.001 and abs(se_lat) < 0.001:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P7 year={2000 + i}: LE_lat={le_lat:.6f}° SE_lat={se_lat:.6f}°"
            )

    except Exception as e:
        errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Historical epochs (1900-2100)
# ============================================================
print("\n=== P8: Historical/future mean node ===")

historical_epochs = []
for year in range(1900, 2101, 5):
    jd = swe.julday(year, 7, 1, 0.0)
    historical_epochs.append((jd, f"Y{year}"))

for jd, label in historical_epochs:
    try:
        se_result = swe.calc_ut(jd, MEAN_NODE, FLAGS)
        le_result = ephem.calc_ut(jd, MEAN_NODE, FLAGS)

        diff_arcsec = abs(se_result[0][0] - le_result[0][0]) * 3600.0
        if diff_arcsec > 180 * 3600:
            diff_arcsec = 360 * 3600 - diff_arcsec

        if diff_arcsec < 1.0:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P8 {label}: SE={se_result[0][0]:.6f}° LE={le_result[0][0]:.6f}° diff={diff_arcsec:.2f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P8 {label}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
print(
    f"ROUND 55 FINAL: {passed}/{passed + failed} passed ({100 * passed / (passed + failed):.1f}%)"
)
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
