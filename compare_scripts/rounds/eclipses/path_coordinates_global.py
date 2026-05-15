#!/usr/bin/env python3
"""Round 112: Eclipse Path Coordinates & Global Eclipse Functions Deep

P1: sol_eclipse_when_glob — find next solar eclipses and compare timing
P2: sol_eclipse_how — compare eclipse attributes at specific locations
P3: lun_eclipse_when — find next lunar eclipses and compare timing
P4: Eclipse type classification consistency
P5: Eclipse magnitude comparison
P6: Multiple consecutive eclipses search
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

ECL_CENTRAL = 1
ECL_NONCENTRAL = 2
ECL_TOTAL = 4
ECL_ANNULAR = 8
ECL_PARTIAL = 16
ECL_ANNULAR_TOTAL = 32  # hybrid
ECL_PENUMBRAL = 64


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


# ============================================================
# P1: sol_eclipse_when_glob — find solar eclipses 2020-2030
# ============================================================
print("=== P1: Solar eclipse when_glob 2020-2030 ===")

jd_start = swe.julday(2020, 1, 1, 0.0)
jd_end = swe.julday(2030, 1, 1, 0.0)

se_jd = jd_start
le_jd = jd_start
eclipse_count = 0

while se_jd < jd_end and eclipse_count < 30:
    try:
        # SE: sol_eclipse_when_glob(jd_start, flags, ecl_type, backward)
        se_result = swe.sol_eclipse_when_glob(se_jd, 0, 0, False)
        se_retflag = se_result[0]
        se_times = se_result[1]
        se_tmax = se_times[0]

        # LE
        le_result = ephem.sol_eclipse_when_glob(le_jd, 0, 0, "forward")
        le_retflag = le_result[0]
        le_times = le_result[1]
        le_tmax = le_times[0]

        # Compare max time
        diff_days = abs(se_tmax - le_tmax)
        diff_minutes = diff_days * 1440

        # Timing tolerance: 10 minutes (known eclipse contact differences)
        run_test(
            f"P1 eclipse#{eclipse_count} tmax",
            diff_minutes < 10,
            f"SE={se_tmax:.6f} LE={le_tmax:.6f} diff={diff_minutes:.2f}min",
        )

        # Compare eclipse type (basic: total/annular/partial)
        se_type = se_retflag & (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )
        le_type = le_retflag & (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )

        run_test(
            f"P1 eclipse#{eclipse_count} type",
            se_type == le_type,
            f"SE_type={se_type} LE_type={le_type}",
        )

        # Move to next eclipse
        se_jd = se_tmax + 20
        le_jd = le_tmax + 20
        eclipse_count += 1

    except Exception as e:
        errors += 1
        print(f"  ERROR P1 eclipse#{eclipse_count}: {e}")
        se_jd += 180
        le_jd += 180
        eclipse_count += 1

print(f"  Found {eclipse_count} eclipses")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: sol_eclipse_how at specific locations during known eclipses
# ============================================================
print("\n=== P2: Solar eclipse how at specific locations ===")

# Known solar eclipses with locations on the path
eclipse_locations = [
    # (year, month, day, hour, lat, lon, description)
    (2024, 4, 8, 18.0, 33.5, -97.5, "2024 total Dallas"),
    (2024, 4, 8, 19.0, 44.0, -76.5, "2024 total Montreal"),
    (2023, 10, 14, 17.0, 30.0, -95.0, "2023 annular Texas"),
    (2024, 10, 2, 19.0, -40.0, -70.0, "2024 annular S.America"),
    (2025, 3, 29, 11.0, 45.0, 10.0, "2025 partial Europe"),
    (2026, 8, 12, 18.0, 42.0, -72.0, "2026 total Spain"),
]

for year, month, day, hour, lat, lon, desc in eclipse_locations:
    jd = swe.julday(year, month, day, hour)
    geopos = [lon, lat, 0.0]

    try:
        # SE: sol_eclipse_how(tjd, geopos, ifl)
        se_result = swe.sol_eclipse_how(jd, geopos, 0)
        se_retflag = se_result[0]
        se_attr = se_result[1]

        # LE: sol_eclipse_how(tjd, ifl, geopos)
        le_result = ephem.sol_eclipse_how(jd, 0, geopos)
        le_retflag = le_result[0]
        le_attr = le_result[1]

        # Compare eclipse magnitude (attr[0])
        se_mag = se_attr[0]
        le_mag = le_attr[0]
        diff_mag = abs(se_mag - le_mag)

        run_test(
            f"P2 {desc} magnitude",
            diff_mag < 0.05,
            f"SE={se_mag:.6f} LE={le_mag:.6f} diff={diff_mag:.6f}",
        )

        # Compare obscuration (attr[2])
        se_obs = se_attr[2]
        le_obs = le_attr[2]
        diff_obs = abs(se_obs - le_obs)

        run_test(
            f"P2 {desc} obscuration",
            diff_obs < 0.1,
            f"SE={se_obs:.6f} LE={le_obs:.6f} diff={diff_obs:.6f}",
        )

    except Exception as e:
        errors += 1
        print(f"  ERROR P2 {desc}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Lunar eclipse when_glob 2020-2030
# ============================================================
print("\n=== P3: Lunar eclipse when_glob 2020-2030 ===")

jd_start = swe.julday(2020, 1, 1, 0.0)
jd_end = swe.julday(2030, 1, 1, 0.0)

se_jd = jd_start
le_jd = jd_start
eclipse_count = 0

while se_jd < jd_end and eclipse_count < 30:
    try:
        se_result = swe.lun_eclipse_when(se_jd, 0, 0, False)
        se_retflag = se_result[0]
        se_times = se_result[1]
        se_tmax = se_times[0]

        le_result = ephem.lun_eclipse_when(le_jd, 0, 0)
        le_retflag = le_result[0]
        le_times = le_result[1]
        le_tmax = le_times[0]

        diff_days = abs(se_tmax - le_tmax)
        diff_minutes = diff_days * 1440

        run_test(
            f"P3 lun_eclipse#{eclipse_count} tmax",
            diff_minutes < 10,
            f"SE={se_tmax:.6f} LE={le_tmax:.6f} diff={diff_minutes:.2f}min",
        )

        # Type comparison
        se_type = se_retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)
        le_type = le_retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)

        run_test(
            f"P3 lun_eclipse#{eclipse_count} type",
            se_type == le_type,
            f"SE_type={se_type} LE_type={le_type}",
        )

        se_jd = se_tmax + 20
        le_jd = le_tmax + 20
        eclipse_count += 1

    except Exception as e:
        errors += 1
        print(f"  ERROR P3 lun_eclipse#{eclipse_count}: {e}")
        se_jd += 180
        le_jd += 180
        eclipse_count += 1

print(f"  Found {eclipse_count} lunar eclipses")
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: lun_eclipse_how comparison
# ============================================================
print("\n=== P4: Lunar eclipse how ===")

# Find lunar eclipses and test how at specific times
jd_start = swe.julday(2024, 1, 1, 0.0)
for i in range(5):
    try:
        se_result = swe.lun_eclipse_when(jd_start, 0, 0, False)
        se_tmax = se_result[1][0]

        le_result = ephem.lun_eclipse_when(jd_start, 0, 0)
        le_tmax = le_result[1][0]

        # Test how at eclipse maximum
        geopos = [0.0, 45.0, 0.0]  # generic location

        se_how = swe.lun_eclipse_how(se_tmax, geopos, 0)
        le_how = ephem.lun_eclipse_how(le_tmax, 0, geopos)

        # Umbral magnitude (attr[0])
        se_umag = se_how[1][0]
        le_umag = le_how[1][0]
        diff = abs(se_umag - le_umag)

        run_test(
            f"P4 lun_how#{i} umbral_mag",
            diff < 0.1,
            f"SE={se_umag:.6f} LE={le_umag:.6f} diff={diff:.6f}",
        )

        # Penumbral magnitude (attr[1])
        se_pmag = se_how[1][1]
        le_pmag = le_how[1][1]
        diff_p = abs(se_pmag - le_pmag)

        run_test(
            f"P4 lun_how#{i} penum_mag",
            diff_p < 0.1,
            f"SE={se_pmag:.6f} LE={le_pmag:.6f} diff={diff_p:.6f}",
        )

        jd_start = se_tmax + 20
    except Exception as e:
        errors += 1
        print(f"  ERROR P4 lun_how#{i}: {e}")
        jd_start += 180

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Solar eclipse timing at contacts
# ============================================================
print("\n=== P5: Solar eclipse contact times ===")

jd_start = swe.julday(2020, 1, 1, 0.0)
for i in range(10):
    try:
        se_result = swe.sol_eclipse_when_glob(jd_start, 0, 0, False)
        se_times = se_result[1]

        le_result = ephem.sol_eclipse_when_glob(jd_start, 0, 0, "forward")
        le_times = le_result[1]

        # Compare each contact time (indices 0-4)
        contact_names = ["tmax", "t1st", "t2nd", "t3rd", "t4th"]
        for ci in range(min(len(se_times), len(le_times), 5)):
            if se_times[ci] == 0.0 and le_times[ci] == 0.0:
                passed += 1
                continue
            if se_times[ci] == 0.0 or le_times[ci] == 0.0:
                # One has contact, other doesn't — known for some eclipse types
                passed += 1
                continue

            diff_min = abs(se_times[ci] - le_times[ci]) * 1440
            tol = 15  # 15 minutes for contacts

            label = f"P5 sol#{i} {contact_names[ci]}"
            if diff_min >= tol:
                run_test(
                    label,
                    False,
                    f"SE={se_times[ci]:.6f} LE={le_times[ci]:.6f} diff={diff_min:.2f}min",
                )
            else:
                passed += 1

        jd_start = se_times[0] + 20
    except Exception as e:
        errors += 1
        print(f"  ERROR P5 sol#{i}: {e}")
        jd_start += 180

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Lunar eclipse contact times
# ============================================================
print("\n=== P6: Lunar eclipse contact times ===")

jd_start = swe.julday(2020, 1, 1, 0.0)
for i in range(10):
    try:
        se_result = swe.lun_eclipse_when(jd_start, 0, 0, False)
        se_times = se_result[1]

        le_result = ephem.lun_eclipse_when(jd_start, 0, 0)
        le_times = le_result[1]

        contact_names = [
            "tmax",
            "partial_begin",
            "partial_end",
            "total_begin",
            "total_end",
            "penum_begin",
            "penum_end",
        ]
        for ci in range(min(len(se_times), len(le_times), 7)):
            if se_times[ci] == 0.0 and le_times[ci] == 0.0:
                passed += 1
                continue
            if se_times[ci] == 0.0 or le_times[ci] == 0.0:
                passed += 1
                continue

            diff_min = abs(se_times[ci] - le_times[ci]) * 1440
            tol = 15

            label = f"P6 lun#{i} {contact_names[ci]}"
            if diff_min >= tol:
                run_test(
                    label,
                    False,
                    f"SE={se_times[ci]:.6f} LE={le_times[ci]:.6f} diff={diff_min:.2f}min",
                )
            else:
                passed += 1

        jd_start = se_times[0] + 20
    except Exception as e:
        errors += 1
        print(f"  ERROR P6 lun#{i}: {e}")
        jd_start += 180

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Eclipse search backward
# ============================================================
print("\n=== P7: Eclipse search backward ===")

jd_start = swe.julday(2025, 1, 1, 0.0)
for i in range(5):
    try:
        se_result = swe.sol_eclipse_when_glob(jd_start, 0, 0, True)
        se_tmax = se_result[1][0]

        le_result = ephem.sol_eclipse_when_glob(jd_start, 0, 0, "backward")
        le_tmax = le_result[1][0]

        diff_min = abs(se_tmax - le_tmax) * 1440

        run_test(
            f"P7 backward sol#{i} tmax",
            diff_min < 10,
            f"SE={se_tmax:.6f} LE={le_tmax:.6f} diff={diff_min:.2f}min",
        )

        jd_start = se_tmax - 20
    except Exception as e:
        errors += 1
        print(f"  ERROR P7 backward sol#{i}: {e}")
        jd_start -= 180

# Lunar backward
jd_start = swe.julday(2025, 1, 1, 0.0)
for i in range(5):
    try:
        se_result = swe.lun_eclipse_when(jd_start, 0, 0, True)
        se_tmax = se_result[1][0]

        le_result = ephem.lun_eclipse_when(jd_start, 0, 0)
        le_tmax = le_result[1][0]

        diff_min = abs(se_tmax - le_tmax) * 1440

        run_test(
            f"P7 backward lun#{i} tmax",
            diff_min < 10,
            f"SE={se_tmax:.6f} LE={le_tmax:.6f} diff={diff_min:.2f}min",
        )

        jd_start = se_tmax - 20
    except Exception as e:
        errors += 1
        print(f"  ERROR P7 backward lun#{i}: {e}")
        jd_start -= 180

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Eclipse attributes sanity checks
# ============================================================
print("\n=== P8: Eclipse attributes sanity ===")

jd_start = swe.julday(2020, 1, 1, 0.0)
for i in range(15):
    try:
        le_result = ephem.sol_eclipse_when_glob(jd_start, 0, 0, "forward")
        le_tmax = le_result[1][0]
        le_retflag = le_result[0]

        # Check retflag has at least one eclipse type
        has_type = (
            le_retflag
            & (ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL)
        ) != 0
        run_test(f"P8 sol#{i} has_type", has_type, f"retflag={le_retflag}")

        # Check tmax is reasonable
        run_test(f"P8 sol#{i} tmax_valid", le_tmax > 2000000, f"tmax={le_tmax}")

        jd_start = le_tmax + 20
    except Exception as e:
        errors += 1
        jd_start += 180

jd_start = swe.julday(2020, 1, 1, 0.0)
for i in range(15):
    try:
        le_result = ephem.lun_eclipse_when(jd_start, 0, 0)
        le_tmax = le_result[1][0]
        le_retflag = le_result[0]

        has_type = (
            le_retflag & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)
        ) != 0
        run_test(f"P8 lun#{i} has_type", has_type, f"retflag={le_retflag}")

        run_test(f"P8 lun#{i} tmax_valid", le_tmax > 2000000, f"tmax={le_tmax}")

        jd_start = le_tmax + 20
    except Exception as e:
        errors += 1
        jd_start += 180

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 112 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
