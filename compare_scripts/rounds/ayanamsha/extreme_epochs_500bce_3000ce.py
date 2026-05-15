#!/usr/bin/env python3
"""Round 82: Ayanamsha Consistency at Extreme Epochs

Tests all ayanamsha modes across a wide range of dates (500 BCE to 3000 CE)
to verify consistency and detect any discontinuities or overflow issues.
Also tests that sidereal planet positions are consistent with the ayanamsha offset.
"""

from __future__ import annotations

import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

print("=" * 70)
print("ROUND 82: Ayanamsha Consistency at Extreme Epochs")
print("=" * 70)

# Ayanamsha modes to test (most common ones)
ayanamsha_modes = [
    (0, "FAGAN_BRADLEY"),
    (1, "LAHIRI"),
    (2, "DELUCE"),
    (3, "RAMAN"),
    (4, "USHASHASHI"),
    (5, "KRISHNAMURTI"),
    (6, "DJWHAL_KHUL"),
    (7, "YUKTESHWAR"),
    (8, "JN_BHASIN"),
    (9, "BABYL_KUGLER1"),
    (10, "BABYL_KUGLER2"),
    (11, "BABYL_KUGLER3"),
    (12, "BABYL_HUBER"),
    (13, "BABYL_ETPSC"),
    (14, "ALDEBARAN_15TAU"),
    (15, "HIPPARCHOS"),
    (16, "SASSANIAN"),
    (17, "GALCENT_0SAG"),
    (21, "GALCENT_COCHRANE"),
    (22, "GALEQU_IAU1958"),
    (23, "GALEQU_TRUE"),
    (27, "TRUE_CITRA"),
    (28, "TRUE_REVATI"),
    (29, "TRUE_PUSHYA"),
    (30, "GALCENT_RGILBRAND"),
    (33, "GALEQU_MULA"),
    (35, "GALCENT_MULA_WILHELM"),
    (36, "ARYABHATA_522"),
    (39, "TRUE_SHEORAN"),
    (40, "GALCENT_COCHRANE_2"),
    (42, "VALENS_MOON"),
]

# ============================================================
# P1: Ayanamsha values across 2000 years
# ============================================================
print("\n=== P1: Ayanamsha values SE vs LE (1000 BCE - 3000 CE) ===")

# Test dates spanning wide range
test_jds = []
for year in [
    -500,
    0,
    500,
    1000,
    1500,
    1600,
    1700,
    1800,
    1900,
    1950,
    2000,
    2025,
    2050,
    2100,
    2200,
    2500,
]:
    jd = swe.julday(year, 1, 1, 12.0)
    test_jds.append((year, jd))

for mode_id, mode_name in ayanamsha_modes:
    mode_pass = 0
    mode_fail = 0
    swe.set_sid_mode(mode_id)
    ephem.set_sid_mode(mode_id, 0.0, 0.0)

    for year, jd in test_jds:
        label = f"{mode_name}(#{mode_id}) y={year}"
        try:
            se_aya = swe.get_ayanamsa_ut(jd)
            le_aya = ephem.get_ayanamsa_ut(jd)
            diff_arcsec = abs(se_aya - le_aya) * 3600.0
            # Known ~14" systematic offset for many modes
            if diff_arcsec < 20.0:
                mode_pass += 1
                passed += 1
            else:
                mode_fail += 1
                failed += 1
                if mode_fail <= 2:
                    print(
                        f'  FAIL {label}: SE={se_aya:.6f} LE={le_aya:.6f} diff={diff_arcsec:.1f}"'
                    )
        except Exception as e:
            errors += 1
            if "range" not in str(e).lower():
                print(f"  ERROR {label}: {e}")

    if mode_fail > 2:
        print(f"  ... {mode_name}: {mode_fail} total failures")

# Reset to default
swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0.0, 0.0)

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Ayanamsha monotonicity (should increase over time for precession)
# ============================================================
print("\n=== P2: Ayanamsha monotonicity check ===")

for mode_id, mode_name in [(0, "FAGAN_BRADLEY"), (1, "LAHIRI"), (27, "TRUE_CITRA")]:
    ephem.set_sid_mode(mode_id, 0.0, 0.0)
    prev_aya = None
    mono_pass = 0
    mono_fail = 0
    for year in range(1900, 2101, 10):
        jd = swe.julday(year, 1, 1, 12.0)
        le_aya = ephem.get_ayanamsa_ut(jd)
        if prev_aya is not None:
            # Precession makes ayanamsha increase by ~50.3"/year
            diff = le_aya - prev_aya
            if diff > 0:  # Should be positive (increasing)
                mono_pass += 1
            else:
                mono_fail += 1
                print(
                    f"  FAIL {mode_name} monotonicity at {year}: aya={le_aya:.6f} prev={prev_aya:.6f}"
                )
        prev_aya = le_aya

    passed += mono_pass
    failed += mono_fail
    print(f"  {mode_name}: {mono_pass} passed, {mono_fail} failed")

ephem.set_sid_mode(0, 0.0, 0.0)
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Sidereal positions = tropical - ayanamsha
# ============================================================
print("\n=== P3: Sidereal = tropical - ayanamsha consistency ===")

FLG_SWIEPH = 2
FLG_SPEED = 256
FLG_SIDEREAL = 65536

for mode_id, mode_name in [(0, "FAGAN_BRADLEY"), (1, "LAHIRI"), (27, "TRUE_CITRA")]:
    for year in [1950, 2000, 2025]:
        jd = swe.julday(year, 7, 1, 12.0)
        for body, body_name in [
            (0, "Sun"),
            (1, "Moon"),
            (2, "Mercury"),
            (4, "Mars"),
            (5, "Jupiter"),
        ]:
            label = f"{mode_name} {body_name} {year}"
            try:
                # Get tropical position
                ephem.set_sid_mode(0, 0.0, 0.0)
                le_trop = ephem.calc_ut(jd, body, FLG_SWIEPH | FLG_SPEED)

                # Get sidereal position
                ephem.set_sid_mode(mode_id, 0.0, 0.0)
                le_sid = ephem.calc_ut(
                    jd, body, FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL
                )

                # Get ayanamsha
                le_aya = ephem.get_ayanamsa_ut(jd)

                # Sidereal longitude should be tropical - ayanamsha
                expected_sid = (le_trop[0] - le_aya) % 360.0
                actual_sid = le_sid[0]
                diff = abs(expected_sid - actual_sid) * 3600.0
                if diff > 180 * 3600:
                    diff = 360 * 3600 - diff

                if diff < 1.0:  # 1" tolerance
                    passed += 1
                else:
                    failed += 1
                    print(
                        f'  FAIL {label}: expected={expected_sid:.6f} actual={actual_sid:.6f} diff={diff:.2f}"'
                    )

            except Exception as e:
                errors += 1
                print(f"  ERROR {label}: {e}")

ephem.set_sid_mode(0, 0.0, 0.0)
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Ayanamsha rate of change (~50.3"/year)
# ============================================================
print("\n=== P4: Ayanamsha precession rate ===")

ephem.set_sid_mode(1, 0.0, 0.0)  # Lahiri
swe.set_sid_mode(1)

for year in [1900, 1950, 2000, 2050, 2100]:
    jd1 = swe.julday(year, 1, 1, 12.0)
    jd2 = swe.julday(year + 1, 1, 1, 12.0)

    le_aya1 = ephem.get_ayanamsa_ut(jd1)
    le_aya2 = ephem.get_ayanamsa_ut(jd2)

    rate_arcsec = (le_aya2 - le_aya1) * 3600.0
    label = f"precession_rate {year}"

    # IAU precession rate is ~50.29"/yr, varies slightly
    if 49.5 < rate_arcsec < 51.5:
        passed += 1
    else:
        failed += 1
        print(f'  FAIL {label}: rate={rate_arcsec:.3f}"/yr (expected ~50.3)')

    # Also compare SE rate
    se_aya1 = swe.get_ayanamsa_ut(jd1)
    se_aya2 = swe.get_ayanamsa_ut(jd2)
    se_rate = (se_aya2 - se_aya1) * 3600.0

    if abs(rate_arcsec - se_rate) < 1.0:  # rates should match within 1"/yr
        passed += 1
    else:
        failed += 1
        print(
            f'  FAIL {label} rate_diff: LE={rate_arcsec:.3f} SE={se_rate:.3f} diff={abs(rate_arcsec - se_rate):.3f}"/yr'
        )

ephem.set_sid_mode(0, 0.0, 0.0)
swe.set_sid_mode(0)
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Custom ayanamsha (user-defined t0/ayan_t0)
# ============================================================
print("\n=== P5: Custom ayanamsha mode ===")

# SIDM_USER = 255
custom_t0 = 2451545.0  # J2000
custom_ayan = 23.5  # Custom initial value

swe.set_sid_mode(255, custom_t0, custom_ayan)
ephem.set_sid_mode(255, custom_t0, custom_ayan)

for year in [1950, 1975, 2000, 2025, 2050]:
    jd = swe.julday(year, 1, 1, 12.0)
    label = f"custom_aya y={year}"
    try:
        se_aya = swe.get_ayanamsa_ut(jd)
        le_aya = ephem.get_ayanamsa_ut(jd)
        diff_arcsec = abs(se_aya - le_aya) * 3600.0
        if diff_arcsec < 20.0:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL {label}: SE={se_aya:.6f} LE={le_aya:.6f} diff={diff_arcsec:.1f}"'
            )
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

# At t0, ayanamsha should equal custom_ayan
le_at_t0 = ephem.get_ayanamsa_ut(custom_t0)
if abs(le_at_t0 - custom_ayan) < 0.01:
    passed += 1
else:
    failed += 1
    print(f"  FAIL custom at t0: expected {custom_ayan} got {le_at_t0:.6f}")

swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0.0, 0.0)
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: True ayanamsha modes (star-based) stability
# ============================================================
print("\n=== P6: True ayanamsha (star-based) stability ===")

true_modes = [
    (27, "TRUE_CITRA"),
    (28, "TRUE_REVATI"),
    (29, "TRUE_PUSHYA"),
    (39, "TRUE_SHEORAN"),
]

for mode_id, mode_name in true_modes:
    swe.set_sid_mode(mode_id)
    ephem.set_sid_mode(mode_id, 0.0, 0.0)

    for year in [1900, 1950, 2000, 2020, 2025]:
        jd = swe.julday(year, 1, 1, 12.0)
        label = f"{mode_name} y={year}"
        try:
            se_aya = swe.get_ayanamsa_ut(jd)
            le_aya = ephem.get_ayanamsa_ut(jd)
            diff_arcsec = abs(se_aya - le_aya) * 3600.0
            if diff_arcsec < 20.0:
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: SE={se_aya:.6f} LE={le_aya:.6f} diff={diff_arcsec:.1f}"'
                )
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0.0, 0.0)
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Ayanamsha SE vs LE comparison at J2000 (reference point)
# ============================================================
print("\n=== P7: Ayanamsha at J2000 reference point ===")

jd_j2000 = 2451545.0
for mode_id, mode_name in ayanamsha_modes:
    swe.set_sid_mode(mode_id)
    ephem.set_sid_mode(mode_id, 0.0, 0.0)
    label = f"{mode_name}(#{mode_id}) at J2000"
    try:
        se_aya = swe.get_ayanamsa_ut(jd_j2000)
        le_aya = ephem.get_ayanamsa_ut(jd_j2000)
        diff_arcsec = abs(se_aya - le_aya) * 3600.0
        if diff_arcsec < 20.0:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL {label}: SE={se_aya:.6f} LE={le_aya:.6f} diff={diff_arcsec:.1f}"'
            )
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0.0, 0.0)
print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# Summary
# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 82 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
