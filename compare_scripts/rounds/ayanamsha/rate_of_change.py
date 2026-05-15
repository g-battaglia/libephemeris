#!/usr/bin/env python3
"""Round 109: Ayanamsha Rate of Change Deep

Tests ayanamsha values AND their time derivatives across all sidereal modes.
P1: get_ayanamsa_ex_ut for all modes at multiple epochs
P2: Ayanamsha rate via finite difference vs analytical (where available)
P3: Ayanamsha continuity (no jumps across year boundaries)
P4: Ayanamsha monotonicity (should be monotonically increasing for most modes)
P5: get_ayanamsa_ex_ut with different flags (J2000, NONUT)
P6: Sidereal position consistency (planet_sidereal = planet_tropical - ayanamsha)
P7: All 43+ ayanamsha modes comparison at J2000
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
FLG_SIDEREAL = 65536
FLG_J2000 = 32
FLG_NONUT = 64

# All ayanamsha modes (0-42+)
AYANAMSHA_MODES = {
    0: "FAGAN_BRADLEY",
    1: "LAHIRI",
    2: "DELUCE",
    3: "RAMAN",
    4: "USHASHASHI",
    5: "KRISHNAMURTI",
    6: "DJWHAL_KHUL",
    7: "YUKTESHWAR",
    8: "JN_BHASIN",
    9: "BABYL_KUGLER1",
    10: "BABYL_KUGLER2",
    11: "BABYL_KUGLER3",
    12: "BABYL_HUBER",
    13: "BABYL_ETPSC",
    14: "ALDEBARAN_15TAU",
    15: "HIPPARCHOS",
    16: "SASSANIAN",
    17: "GALCENT_0SAG",
    18: "J2000",
    19: "J1900",
    20: "B1950",
    21: "SURYASIDDHANTA",
    22: "SURYASIDDHANTA_MSUN",
    23: "ARYABHATA",
    24: "ARYABHATA_MSUN",
    25: "SS_REVATI",
    26: "SS_CITRA",
    27: "TRUE_CITRA",
    28: "TRUE_REVATI",
    29: "TRUE_PUSHYA",
    30: "GALCENT_RGILBRAND",
    31: "GALEQU_IAU1958",
    32: "GALEQU_TRUE",
    33: "GALEQU_MULA",
    34: "GALCENT_MULA_WILHELM",
    35: "ARYABHATA_522",
    36: "BABYL_BRITTON",
    37: "TRUE_SHEORAN",
    38: "GALCENT_COCHRANE",
    39: "GALCENT_FIORENZA",
    40: "VALENS_MOON",
    41: "LAHIRI_1940",
    42: "LAHIRI_VP285",
    43: "KRISHNAMURTI_VP291",
    44: "LAHIRI_ICRC",
}


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


# ============================================================
# P1: Ayanamsha values at multiple epochs for all modes
# ============================================================
print("=== P1: Ayanamsha values all modes at key epochs ===")

test_epochs = [
    (swe.julday(1900, 1, 1, 0.0), "Y1900"),
    (swe.julday(1950, 1, 1, 0.0), "Y1950"),
    (2451545.0, "J2000"),
    (swe.julday(2000, 1, 1, 12.0), "Y2000"),
    (swe.julday(2024, 1, 1, 0.0), "Y2024"),
    (swe.julday(2050, 1, 1, 0.0), "Y2050"),
]

for mode_id, mode_name in AYANAMSHA_MODES.items():
    for jd, epoch_label in test_epochs:
        try:
            swe.set_sid_mode(mode_id)
            ephem.set_sid_mode(mode_id, 0, 0)

            se_aya = swe.get_ayanamsa_ut(jd)
            le_aya = ephem.get_ayanamsa_ut(jd)

            diff = abs(se_aya - le_aya)
            diff_arcsec = diff * 3600

            # Known ~14" systematic offset for some modes
            # Use 20" tolerance to accommodate known sidereal model differences
            tol_arcsec = 20.0

            # Galactic center modes at extreme dates can be much worse
            if mode_id in (17, 30, 35, 38) and epoch_label in ("Y1900",):
                tol_arcsec = 800.0

            label = f"P1 {mode_name} {epoch_label}"
            run_test(
                label,
                diff_arcsec < tol_arcsec,
                f'SE={se_aya:.8f}° LE={le_aya:.8f}° diff={diff_arcsec:.2f}"',
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P1 {mode_name} {epoch_label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Ayanamsha rate via finite difference
# Compare LE rate at each epoch using central differences
# ============================================================
print("\n=== P2: Ayanamsha rate of change (finite difference) ===")

DT = 1.0  # 1 day for finite difference
rate_epochs = [
    (2451545.0, "J2000"),
    (swe.julday(2024, 6, 15, 12.0), "2024-mid"),
    (swe.julday(1950, 6, 15, 12.0), "1950-mid"),
]

# Test a subset of modes for rate
RATE_MODES = {
    0: "FAGAN_BRADLEY",
    1: "LAHIRI",
    3: "RAMAN",
    5: "KRISHNAMURTI",
    14: "ALDEBARAN_15TAU",
    18: "J2000",
    27: "TRUE_CITRA",
    28: "TRUE_REVATI",
    29: "TRUE_PUSHYA",
    41: "LAHIRI_1940",
}

for mode_id, mode_name in RATE_MODES.items():
    for jd, epoch_label in rate_epochs:
        try:
            swe.set_sid_mode(mode_id)
            ephem.set_sid_mode(mode_id, 0, 0)

            # SE rate via finite diff
            se_minus = swe.get_ayanamsa_ut(jd - DT)
            se_plus = swe.get_ayanamsa_ut(jd + DT)
            se_rate = (se_plus - se_minus) / (2 * DT)  # deg/day

            # LE rate via finite diff
            le_minus = ephem.get_ayanamsa_ut(jd - DT)
            le_plus = ephem.get_ayanamsa_ut(jd + DT)
            le_rate = (le_plus - le_minus) / (2 * DT)  # deg/day

            # Convert to arcsec/year for readability
            se_rate_arcsec_yr = se_rate * 3600 * 365.25
            le_rate_arcsec_yr = le_rate * 3600 * 365.25

            diff_arcsec_yr = abs(se_rate_arcsec_yr - le_rate_arcsec_yr)

            # Rate tolerance: 0.5 arcsec/year (very generous — precession rate ~50"/yr)
            tol = 0.5

            # TRUE_* modes have nutation-dependent rates, allow wider
            if mode_id in (27, 28, 29, 32, 33, 37, 40):
                tol = 2.0

            label = f"P2 {mode_name} {epoch_label}"
            run_test(
                label,
                diff_arcsec_yr < tol,
                f'SE_rate={se_rate_arcsec_yr:.4f}"/yr LE_rate={le_rate_arcsec_yr:.4f}"/yr diff={diff_arcsec_yr:.4f}"/yr',
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P2 {mode_name} {epoch_label}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Ayanamsha continuity across year boundaries
# Check no discontinuities
# ============================================================
print("\n=== P3: Ayanamsha continuity across year boundaries ===")

CONT_MODES = {0: "FAGAN_BRADLEY", 1: "LAHIRI", 27: "TRUE_CITRA"}

for mode_id, mode_name in CONT_MODES.items():
    ephem.set_sid_mode(mode_id, 0, 0)

    for year in range(1990, 2031):
        jd_end = swe.julday(year, 12, 31, 23.0)
        jd_start = swe.julday(year + 1, 1, 1, 1.0)
        try:
            le_end = ephem.get_ayanamsa_ut(jd_end)
            le_start = ephem.get_ayanamsa_ut(jd_start)

            # 2 hours apart, rate ~50"/yr → expected change ~0.011"
            expected_max_change = 0.001  # degrees (~3.6")
            actual_change = abs(le_start - le_end)

            label = f"P3 {mode_name} {year}/{year + 1}"
            run_test(
                label,
                actual_change < expected_max_change,
                f'end={le_end:.8f} start={le_start:.8f} jump={actual_change * 3600:.2f}"',
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P3 {mode_name} {year}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Ayanamsha monotonicity (most modes should be monotonically increasing)
# ============================================================
print("\n=== P4: Ayanamsha monotonicity ===")

# Non-monotonic modes: TRUE_* modes may oscillate due to nutation
# Most mean-precession modes should be strictly increasing
MONO_MODES = {
    0: "FAGAN_BRADLEY",
    1: "LAHIRI",
    3: "RAMAN",
    5: "KRISHNAMURTI",
    14: "ALDEBARAN_15TAU",
    18: "J2000",
}

for mode_id, mode_name in MONO_MODES.items():
    ephem.set_sid_mode(mode_id, 0, 0)

    prev_aya = None
    monotonic = True
    violation_detail = ""

    for year in range(1900, 2101, 5):
        jd = swe.julday(year, 1, 1, 12.0)
        try:
            le_aya = ephem.get_ayanamsa_ut(jd)
            if prev_aya is not None:
                if le_aya <= prev_aya:
                    monotonic = False
                    violation_detail = f"Y{year}: {le_aya:.8f} <= prev {prev_aya:.8f}"
                    break
            prev_aya = le_aya
        except Exception as e:
            errors += 1

    label = f"P4 {mode_name} monotonic"
    run_test(label, monotonic, violation_detail)

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Ayanamsha with different flags (get_ayanamsa_ex_ut)
# ============================================================
print("\n=== P5: Ayanamsha with flags (get_ayanamsa_ex_ut) ===")

FLAG_COMBOS = [
    (0, "default"),
    (FLG_J2000 | FLG_NONUT, "J2000+NONUT"),
    (FLG_NONUT, "NONUT"),
]

ex_modes = {0: "FAGAN_BRADLEY", 1: "LAHIRI", 27: "TRUE_CITRA", 18: "J2000"}

for mode_id, mode_name in ex_modes.items():
    for flags, flag_label in FLAG_COMBOS:
        for jd, epoch_label in test_epochs:
            try:
                swe.set_sid_mode(mode_id)
                ephem.set_sid_mode(mode_id, 0, 0)

                # SE: swe.get_ayanamsa_ex_ut(jd, flags) returns (retflag, aya)
                se_result = swe.get_ayanamsa_ex_ut(jd, flags)
                if isinstance(se_result, tuple):
                    se_aya = se_result[1] if len(se_result) >= 2 else se_result[0]
                else:
                    se_aya = se_result

                # LE
                le_result = ephem.get_ayanamsa_ex_ut(jd, flags)
                if isinstance(le_result, tuple):
                    le_aya = le_result[1] if len(le_result) >= 2 else le_result[0]
                else:
                    le_aya = le_result

                diff_arcsec = abs(se_aya - le_aya) * 3600

                tol = 20.0  # 20" known sidereal offset
                if mode_id in (27, 28, 29) and "NONUT" in flag_label:
                    tol = 25.0  # TRUE_* with NONUT may differ more

                label = f"P5 {mode_name} {flag_label} {epoch_label}"
                run_test(
                    label,
                    diff_arcsec < tol,
                    f'SE={se_aya:.8f}° LE={le_aya:.8f}° diff={diff_arcsec:.2f}"',
                )
            except Exception as e:
                errors += 1
                print(f"  ERROR P5 {mode_name} {flag_label} {epoch_label}: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Sidereal position consistency
# planet_sidereal ≈ planet_tropical - ayanamsha (for ecliptic longitude)
# ============================================================
print("\n=== P6: Sidereal position = tropical - ayanamsha consistency ===")

CONSIST_MODES = {0: "FAGAN_BRADLEY", 1: "LAHIRI", 27: "TRUE_CITRA"}
CONSIST_BODIES = {0: "Sun", 1: "Moon", 2: "Mercury", 4: "Mars", 5: "Jupiter"}

for mode_id, mode_name in CONSIST_MODES.items():
    for jd, epoch_label in test_epochs[:4]:
        for body_id, body_name in CONSIST_BODIES.items():
            try:
                ephem.set_sid_mode(mode_id, 0, 0)

                # Tropical position
                le_trop = ephem.calc_ut(jd, body_id, FLG_SPEED)
                trop_lon = le_trop[0][0]

                # Sidereal position
                le_sid = ephem.calc_ut(jd, body_id, FLG_SPEED | FLG_SIDEREAL)
                sid_lon = le_sid[0][0]

                # Ayanamsha
                le_aya = ephem.get_ayanamsa_ut(jd)

                # Compute expected sidereal
                expected_sid = (trop_lon - le_aya) % 360.0

                # Compare
                diff = abs(sid_lon - expected_sid)
                if diff > 180:
                    diff = 360 - diff
                diff_arcsec = diff * 3600

                # Should be very close (same computation internally)
                tol = 1.0  # 1 arcsec tolerance

                label = f"P6 {mode_name} {body_name} {epoch_label}"
                run_test(
                    label,
                    diff_arcsec < tol,
                    f'sid={sid_lon:.8f} expected={expected_sid:.8f} diff={diff_arcsec:.4f}"',
                )
            except Exception as e:
                errors += 1
                print(f"  ERROR P6 {mode_name} {body_name} {epoch_label}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: All modes at J2000 — comparison with SE
# ============================================================
print("\n=== P7: All ayanamsha modes at J2000 ===")

jd_j2000 = 2451545.0
for mode_id, mode_name in AYANAMSHA_MODES.items():
    try:
        swe.set_sid_mode(mode_id)
        ephem.set_sid_mode(mode_id, 0, 0)

        se_aya = swe.get_ayanamsa_ut(jd_j2000)
        le_aya = ephem.get_ayanamsa_ut(jd_j2000)

        diff_arcsec = abs(se_aya - le_aya) * 3600

        # At J2000, most modes should agree well
        tol = 15.0
        if mode_id in (17, 30, 35, 38):
            tol = 500.0  # galactic center modes

        label = f"P7 {mode_name} J2000"
        run_test(
            label,
            diff_arcsec < tol,
            f'SE={se_aya:.8f}° LE={le_aya:.8f}° diff={diff_arcsec:.4f}"',
        )
    except Exception as e:
        errors += 1
        print(f"  ERROR P7 {mode_name}: {e}")

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Ayanamsha rate smoothness (no glitches in derivative)
# ============================================================
print("\n=== P8: Ayanamsha rate smoothness ===")

SMOOTH_MODES = {0: "FAGAN_BRADLEY", 1: "LAHIRI", 14: "ALDEBARAN_15TAU"}

for mode_id, mode_name in SMOOTH_MODES.items():
    ephem.set_sid_mode(mode_id, 0, 0)
    swe.set_sid_mode(mode_id)

    prev_rate = None
    jd_start = swe.julday(2000, 1, 1, 0.0)

    for i in range(120):  # 10 years at monthly intervals
        jd = jd_start + i * 30.4375  # ~monthly
        try:
            # LE rate via finite diff
            dt = 5.0
            le_m = ephem.get_ayanamsa_ut(jd - dt)
            le_p = ephem.get_ayanamsa_ut(jd + dt)
            le_rate = (le_p - le_m) / (2 * dt) * 3600 * 365.25  # arcsec/yr

            # SE rate
            se_m = swe.get_ayanamsa_ut(jd - dt)
            se_p = swe.get_ayanamsa_ut(jd + dt)
            se_rate = (se_p - se_m) / (2 * dt) * 3600 * 365.25

            diff_rate = abs(le_rate - se_rate)

            # Rate should agree within 0.3"/yr
            tol = 0.3

            label = f"P8 {mode_name} month={i}"
            if diff_rate >= tol:
                run_test(
                    label,
                    False,
                    f'LE_rate={le_rate:.4f}"/yr SE_rate={se_rate:.4f}"/yr diff={diff_rate:.4f}"/yr',
                )
            else:
                passed += 1

            # Check smoothness: rate shouldn't jump more than 1"/yr between months
            if prev_rate is not None:
                jump = abs(le_rate - prev_rate)
                if jump > 1.0:
                    run_test(
                        f"P8 {mode_name} smooth month={i}",
                        False,
                        f'rate_jump={jump:.4f}"/yr',
                    )
                else:
                    passed += 1

            prev_rate = le_rate
        except Exception as e:
            errors += 1

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P9: Ayanamsha at century boundaries (long-term stability)
# ============================================================
print("\n=== P9: Ayanamsha at century boundaries ===")

CENTURY_MODES = {0: "FAGAN_BRADLEY", 1: "LAHIRI", 18: "J2000"}

for mode_id, mode_name in CENTURY_MODES.items():
    swe.set_sid_mode(mode_id)
    ephem.set_sid_mode(mode_id, 0, 0)

    for year in range(1600, 2601, 100):
        jd = swe.julday(year, 1, 1, 12.0)
        try:
            se_aya = swe.get_ayanamsa_ut(jd)
            le_aya = ephem.get_ayanamsa_ut(jd)

            diff_arcsec = abs(se_aya - le_aya) * 3600

            # Wider tolerance for extreme dates
            tol = 15.0
            if year < 1800 or year > 2200:
                tol = 30.0

            label = f"P9 {mode_name} Y{year}"
            run_test(
                label,
                diff_arcsec < tol,
                f'SE={se_aya:.6f}° LE={le_aya:.6f}° diff={diff_arcsec:.2f}"',
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P9 {mode_name} Y{year}: {e}")

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P10: Lahiri ayanamsha precision (most commonly used)
# Fine-grained daily test over 2024
# ============================================================
print("\n=== P10: Lahiri daily precision 2024 ===")

swe.set_sid_mode(1)
ephem.set_sid_mode(1, 0, 0)

jd_start = swe.julday(2024, 1, 1, 0.0)
for day in range(365):
    jd = jd_start + day
    try:
        se_aya = swe.get_ayanamsa_ut(jd)
        le_aya = ephem.get_ayanamsa_ut(jd)

        diff_arcsec = abs(se_aya - le_aya) * 3600

        tol = 15.0  # Known ~14" systematic offset

        if diff_arcsec >= tol:
            run_test(
                f"P10 Lahiri day={day}",
                False,
                f'SE={se_aya:.8f}° LE={le_aya:.8f}° diff={diff_arcsec:.4f}"',
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
print(f"ROUND 109 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
