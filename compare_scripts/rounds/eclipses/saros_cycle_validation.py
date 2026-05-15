#!/usr/bin/env python3
"""Round 59: Eclipse Saros Cycle Validation

Verify eclipse prediction against the Saros cycle:
- Saros cycle: 6585.32 days (18 years, 11 days, 8 hours)
- After one Saros, similar eclipses recur
- Compare sol_eclipse_when_glob results between LE and SE
- Verify that eclipses found by LE match Saros periodicity
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLAGS = 0
SAROS_DAYS = 6585.3213  # Saros cycle in days

print("=" * 70)
print("ROUND 59: Eclipse Saros Cycle Validation")
print("=" * 70)


# ============================================================
# P1: Solar eclipse timing comparison (2020-2030)
# ============================================================
print("\n=== P1: Solar eclipse timing 2020-2030 ===")

jd = swe.julday(2020, 1, 1, 0.0)
jd_end = swe.julday(2030, 1, 1, 0.0)

le_eclipses = []
se_eclipses = []

# Find all solar eclipses using LE
jd_search = jd
while jd_search < jd_end:
    try:
        le_result = ephem.sol_eclipse_when_glob(jd_search, FLAGS, 0)
        le_tret = le_result[1] if isinstance(le_result, tuple) else le_result
        le_max = le_tret[0]
        if le_max > 0 and le_max < jd_end:
            le_eclipses.append(
                (le_max, le_result[0] if isinstance(le_result, tuple) else 0)
            )
            jd_search = le_max + 20
        else:
            break
    except Exception as e:
        errors += 1
        print(f"  ERROR P1 LE: {e}")
        break

# Find all solar eclipses using SE
jd_search = jd
while jd_search < jd_end:
    try:
        se_result = swe.sol_eclipse_when_glob(jd_search, FLAGS)
        se_tret = se_result[1]
        se_max = se_tret[0]
        if se_max > 0 and se_max < jd_end:
            se_eclipses.append((se_max, se_result[0]))
            jd_search = se_max + 20
        else:
            break
    except Exception as e:
        errors += 1
        print(f"  ERROR P1 SE: {e}")
        break

print(f"  LE found {len(le_eclipses)} solar eclipses, SE found {len(se_eclipses)}")

if len(le_eclipses) == len(se_eclipses):
    passed += 1
    for i, ((le_max, le_type), (se_max, se_type)) in enumerate(
        zip(le_eclipses, se_eclipses)
    ):
        diff_min = abs(le_max - se_max) * 1440  # minutes
        y, m, d, h = swe.revjul(le_max)
        if diff_min < 30.0:  # Within 30 minutes
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P1 eclipse {i} {y}-{m:02d}-{d:02d}: diff={diff_min:.1f}min")
else:
    failed += 1
    print(f"  FAIL P1: count mismatch LE={len(le_eclipses)} SE={len(se_eclipses)}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Saros cycle periodicity check
# ============================================================
print("\n=== P2: Saros cycle periodicity ===")

# For each eclipse found, check if an eclipse exists ~1 Saros later
jd_2040 = swe.julday(2040, 1, 1, 0.0)

for i, (ecl_jd, ecl_type) in enumerate(le_eclipses[:5]):  # First 5
    expected_next = ecl_jd + SAROS_DAYS
    if expected_next > jd_2040:
        continue
    try:
        # Search for eclipse near expected Saros return
        result = ephem.sol_eclipse_when_glob(expected_next - 15, FLAGS, 0)
        tret = result[1] if isinstance(result, tuple) else result
        found_jd = tret[0]

        diff_hours = abs(found_jd - expected_next) * 24
        y1, m1, d1, h1 = swe.revjul(ecl_jd)
        y2, m2, d2, h2 = swe.revjul(found_jd)

        # Eclipse should occur within ~12 hours of predicted Saros return
        if diff_hours < 24.0:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P2 eclipse {y1}-{m1:02d}-{d1:02d} -> {y2}-{m2:02d}-{d2:02d}: diff={diff_hours:.1f}h"
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P2: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Lunar eclipse timing comparison (2020-2030)
# ============================================================
print("\n=== P3: Lunar eclipse timing 2020-2030 ===")

le_lunar = []
se_lunar = []

jd_search = jd
while jd_search < jd_end:
    try:
        le_result = ephem.lun_eclipse_when(jd_search, FLAGS, 0)
        le_tret = le_result[1] if isinstance(le_result, tuple) else le_result
        le_max = le_tret[0]
        if le_max > 0 and le_max < jd_end:
            le_lunar.append(
                (le_max, le_result[0] if isinstance(le_result, tuple) else 0)
            )
            jd_search = le_max + 20
        else:
            break
    except Exception as e:
        errors += 1
        break

jd_search = jd
while jd_search < jd_end:
    try:
        se_result = swe.lun_eclipse_when(jd_search, FLAGS)
        se_tret = se_result[1]
        se_max = se_tret[0]
        if se_max > 0 and se_max < jd_end:
            se_lunar.append((se_max, se_result[0]))
            jd_search = se_max + 20
        else:
            break
    except Exception as e:
        errors += 1
        break

print(f"  LE found {len(le_lunar)} lunar eclipses, SE found {len(se_lunar)}")

if len(le_lunar) == len(se_lunar):
    passed += 1
    for i, ((le_max, le_type), (se_max, se_type)) in enumerate(zip(le_lunar, se_lunar)):
        diff_min = abs(le_max - se_max) * 1440
        y, m, d, h = swe.revjul(le_max)
        if diff_min < 30.0:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P3 lunar {i} {y}-{m:02d}-{d:02d}: diff={diff_min:.1f}min")
else:
    failed += 1
    print(f"  FAIL P3: count mismatch LE={len(le_lunar)} SE={len(se_lunar)}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Lunar Saros cycle
# ============================================================
print("\n=== P4: Lunar Saros cycle periodicity ===")

for i, (ecl_jd, ecl_type) in enumerate(le_lunar[:5]):
    expected_next = ecl_jd + SAROS_DAYS
    if expected_next > jd_2040:
        continue
    try:
        result = ephem.lun_eclipse_when(expected_next - 15, FLAGS, 0)
        tret = result[1] if isinstance(result, tuple) else result
        found_jd = tret[0]

        diff_hours = abs(found_jd - expected_next) * 24
        y1, m1, d1, h1 = swe.revjul(ecl_jd)
        y2, m2, d2, h2 = swe.revjul(found_jd)

        if diff_hours < 24.0:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P4 lunar {y1}-{m1:02d}-{d1:02d} -> {y2}-{m2:02d}-{d2:02d}: diff={diff_hours:.1f}h"
            )

    except Exception as e:
        errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Eclipse type flags comparison
# ============================================================
print("\n=== P5: Eclipse type flags ===")

for i, ((le_max, le_type), (se_max, se_type)) in enumerate(
    zip(le_eclipses, se_eclipses)
):
    y, m, d, h = swe.revjul(le_max)

    # Check basic type bits match
    # ECL_TOTAL = 1, ECL_ANNULAR = 2, ECL_PARTIAL = 4,
    # ECL_ANNULAR_TOTAL/HYBRID = 8
    le_basic = le_type & 0xF  # Lower 4 bits = eclipse type
    se_basic = se_type & 0xF

    if le_basic == se_basic:
        passed += 1
    else:
        failed += 1
        print(
            f"  FAIL P5 solar {i} {y}-{m:02d}-{d:02d}: LE_type=0x{le_type:x} SE_type=0x{se_type:x}"
        )

for i, ((le_max, le_type), (se_max, se_type)) in enumerate(zip(le_lunar, se_lunar)):
    y, m, d, h = swe.revjul(le_max)
    le_basic = le_type & 0xF
    se_basic = se_type & 0xF

    if le_basic == se_basic:
        passed += 1
    else:
        failed += 1
        print(
            f"  FAIL P5 lunar {i} {y}-{m:02d}-{d:02d}: LE_type=0x{le_type:x} SE_type=0x{se_type:x}"
        )

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Eclipse count per decade (sanity check)
# ============================================================
print("\n=== P6: Eclipse frequency sanity check ===")

# Typically 2-5 solar eclipses per year, ~15-25 per decade
solar_count = len(le_eclipses)
lunar_count = len(le_lunar)

# 2020-2030: 10 years
if 15 <= solar_count <= 30:
    passed += 1
else:
    failed += 1
    print(f"  FAIL P6: {solar_count} solar eclipses in 10yr (expected 15-30)")

if 10 <= lunar_count <= 25:
    passed += 1
else:
    failed += 1
    print(f"  FAIL P6: {lunar_count} lunar eclipses in 10yr (expected 10-25)")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Known eclipses verification
# ============================================================
print("\n=== P7: Known historical eclipses ===")

# Known solar eclipses:
known_solar = [
    (2020, 6, 21, "annular"),  # Annular eclipse 2020-06-21
    (2020, 12, 14, "total"),  # Total eclipse 2020-12-14
    (2021, 6, 10, "annular"),  # Annular eclipse 2021-06-10
    (2021, 12, 4, "total"),  # Total eclipse 2021-12-04
    (2023, 4, 20, "hybrid"),  # Hybrid eclipse 2023-04-20
    (2023, 10, 14, "annular"),  # Annular eclipse 2023-10-14
    (2024, 4, 8, "total"),  # Total eclipse 2024-04-08
    (2024, 10, 2, "annular"),  # Annular eclipse 2024-10-02
]

for year, month, day, expected_type in known_solar:
    jd_search = swe.julday(year, month, day - 5, 0.0)
    try:
        result = ephem.sol_eclipse_when_glob(jd_search, FLAGS, 0)
        tret = result[1] if isinstance(result, tuple) else result
        ecl_jd = tret[0]
        ecl_type = result[0] if isinstance(result, tuple) else 0

        y, m, d, h = swe.revjul(ecl_jd)

        # Date should match within ±1 day
        expected_jd = swe.julday(year, month, day, 0.0)
        diff_days = abs(ecl_jd - expected_jd)

        if diff_days < 1.5:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P7 {year}-{month:02d}-{day:02d}: found at {y}-{m:02d}-{d:02d} (diff={diff_days:.1f}d)"
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P7 {year}-{month:02d}-{day:02d}: {e}")

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Solar eclipse maximum time precision (LE vs SE)
# ============================================================
print("\n=== P8: Eclipse tret array comparison ===")

jd_search = swe.julday(2024, 1, 1, 0.0)
for ecl_num in range(3):
    try:
        le_result = ephem.sol_eclipse_when_glob(jd_search, FLAGS, 0)
        se_result = swe.sol_eclipse_when_glob(jd_search, FLAGS)

        le_tret = le_result[1] if isinstance(le_result, tuple) else le_result
        se_tret = se_result[1]

        # Compare tret values
        labels = [
            "t_max",
            "t_local_noon/begin",
            "t_end",
            "t_begin_total",
            "t_end_total",
        ]
        for j in range(min(5, len(le_tret), len(se_tret))):
            if le_tret[j] == 0.0 and se_tret[j] == 0.0:
                passed += 1
                continue
            if le_tret[j] == 0.0 or se_tret[j] == 0.0:
                # One has the value, other doesn't
                failed += 1
                label = labels[j] if j < len(labels) else f"tret[{j}]"
                print(
                    f"  FAIL P8 ecl{ecl_num} {label}: LE={le_tret[j]:.6f} SE={se_tret[j]:.6f}"
                )
                continue
            diff_min = abs(le_tret[j] - se_tret[j]) * 1440
            label = labels[j] if j < len(labels) else f"tret[{j}]"
            if diff_min < 30.0:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL P8 ecl{ecl_num} {label}: diff={diff_min:.1f}min")

        jd_search = le_tret[0] + 20

    except Exception as e:
        errors += 1
        print(f"  ERROR P8 ecl{ecl_num}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
print(f"ROUND 59 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
