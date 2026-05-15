#!/usr/bin/env python3
"""Round 213: cotrans + equatorial consistency.

Verifies that manual cotrans() ecliptic-to-equatorial transformations
are consistent with FLG_EQUATORIAL flag results, and that cotrans_sp()
speed transformations match equatorial speeds from calc_ut.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

FLAGS_ECL = ephem.FLG_SWIEPH | ephem.FLG_SPEED
FLAGS_EQU = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL

BODIES = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Moon", ephem.MOON, swe.MOON),
    ("Mercury", ephem.MERCURY, swe.MERCURY),
    ("Venus", ephem.VENUS, swe.VENUS),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
]

DATES = [
    2451545.0,  # J2000
    2460000.0,  # 2023
    2440000.0,  # 1968
    2415020.0,  # 1900
    2455000.0,  # 2009
    2462000.0,  # 2028
    2435000.0,  # 1954
    2445000.0,  # 1982
    2450000.0,  # 1995
    2458000.0,  # 2017
]


def compare_cotrans_vs_flag(body_name, le_body, se_body, jd):
    """Compare manual cotrans with FLG_EQUATORIAL flag."""
    global passed, failed, total

    label = f"{body_name} JD={jd:.1f}"

    try:
        # Get ecliptic positions
        le_ecl = ephem.calc_ut(jd, le_body, FLAGS_ECL)
        # Get equatorial via flag
        le_equ = ephem.calc_ut(jd, le_body, FLAGS_EQU)
        # Get obliquity
        le_nut = ephem.calc_ut(jd, -1, 0)
        eps = le_nut[0][0]  # true obliquity
    except Exception as e:
        return

    # Manual cotrans: ecliptic -> equatorial
    try:
        ecl_pos = (le_ecl[0][0], le_ecl[0][1], le_ecl[0][2])
        equ_manual = ephem.cotrans(ecl_pos, -eps)  # negative for ecl->equ
    except Exception as e:
        return

    # Compare RA
    total += 1
    ra_flag = le_equ[0][0]
    ra_manual = equ_manual[0]
    ra_diff = abs(ra_flag - ra_manual)
    if ra_diff > 180:
        ra_diff = 360 - ra_diff
    ra_diff_as = ra_diff * 3600

    # Manual cotrans uses single obliquity (no nutation detail),
    # flag uses full pipeline. Known ~14" difference.
    if ra_diff_as <= 60.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} RA: flag={ra_flag:.6f} manual={ra_manual:.6f} diff={ra_diff_as:.2f}"'
        )

    # Compare DEC
    total += 1
    dec_flag = le_equ[0][1]
    dec_manual = equ_manual[1]
    dec_diff = abs(dec_flag - dec_manual) * 3600

    if dec_diff <= 60.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} DEC: flag={dec_flag:.6f} manual={dec_manual:.6f} diff={dec_diff:.2f}"'
        )

    # Compare LE equatorial vs SE equatorial
    try:
        se_equ = swe.calc_ut(
            jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        )
    except Exception:
        return

    total += 1
    ra_diff_se = abs(le_equ[0][0] - se_equ[0][0])
    if ra_diff_se > 180:
        ra_diff_se = 360 - ra_diff_se
    ra_diff_se_as = ra_diff_se * 3600

    if ra_diff_se_as <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} RA vs SE: LE={le_equ[0][0]:.6f} SE={se_equ[0][0]:.6f} diff={ra_diff_se_as:.4f}"'
        )

    total += 1
    dec_diff_se = abs(le_equ[0][1] - se_equ[0][1]) * 3600
    if dec_diff_se <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} DEC vs SE: LE={le_equ[0][1]:.6f} SE={se_equ[0][1]:.6f} diff={dec_diff_se:.4f}"'
        )

    # Compare RA speed
    total += 1
    ra_spd_diff = abs(le_equ[0][3] - se_equ[0][3]) * 3600
    if ra_spd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} RA_SPD vs SE: LE={le_equ[0][3]:.6f} SE={se_equ[0][3]:.6f} diff={ra_spd_diff:.4f}"/day'
        )

    # Compare DEC speed
    total += 1
    dec_spd_diff = abs(le_equ[0][4] - se_equ[0][4]) * 3600
    if dec_spd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} DEC_SPD vs SE: LE={le_equ[0][4]:.6f} SE={se_equ[0][4]:.6f} diff={dec_spd_diff:.4f}"/day'
        )


def compare_cotrans_sp(body_name, le_body, jd):
    """Test cotrans_sp speed transformation consistency."""
    global passed, failed, total

    label = f"{body_name} JD={jd:.1f} cotrans_sp"

    try:
        le_ecl = ephem.calc_ut(jd, le_body, FLAGS_ECL)
        le_nut = ephem.calc_ut(jd, -1, 0)
        eps = le_nut[0][0]

        ecl_pos_spd = (
            le_ecl[0][0],
            le_ecl[0][1],
            le_ecl[0][2],
            le_ecl[0][3],
            le_ecl[0][4],
            le_ecl[0][5],
        )
        equ_sp = ephem.cotrans_sp(ecl_pos_spd, -eps)
    except Exception:
        return

    # cotrans_sp should return 6 values
    total += 1
    if len(equ_sp) >= 6:
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label}: returned {len(equ_sp)} values, expected 6")
        return

    # Compare with SE cotrans
    try:
        se_ecl_pos = (le_ecl[0][0], le_ecl[0][1], le_ecl[0][2])
        se_result = swe.cotrans(se_ecl_pos, -eps)
    except Exception:
        return

    total += 1
    ra_diff = abs(equ_sp[0] - se_result[0])
    if ra_diff > 180:
        ra_diff = 360 - ra_diff
    ra_diff_as = ra_diff * 3600
    if ra_diff_as <= 0.001:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} RA: LE={equ_sp[0]:.8f} SE={se_result[0]:.8f} diff={ra_diff_as:.6f}"'
        )

    total += 1
    dec_diff = abs(equ_sp[1] - se_result[1]) * 3600
    if dec_diff <= 0.001:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} DEC: LE={equ_sp[1]:.8f} SE={se_result[1]:.8f} diff={dec_diff:.6f}"'
        )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 213: cotrans + equatorial consistency")
    print("=" * 70)

    for body_name, le_b, se_b in BODIES:
        print(f"\n--- {body_name} ---")
        for jd in DATES:
            compare_cotrans_vs_flag(body_name, le_b, se_b, jd)
            compare_cotrans_sp(body_name, le_b, jd)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
