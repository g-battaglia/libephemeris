#!/usr/bin/env python3
"""Round 94: Precession Rates at Extreme Dates

Tests precession by comparing J2000 vs ecliptic-of-date positions across
wide date ranges. Verifies the precession transformation is correct by
checking that the difference between J2000 and of-date grows at ~50.3"/yr.
"""

from __future__ import annotations
import os, sys, math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0
F = 2
S = 256
J = 32
N = 64

print("=" * 70)
print("ROUND 94: Precession Rates at Extreme Dates")
print("=" * 70)

# ============================================================
# P1: Precession offset (of-date - J2000) grows at ~50.3"/yr
# ============================================================
print("\n=== P1: Precession rate from Sun position ===")

for year in range(1950, 2051, 5):
    jd = swe.julday(year, 1, 1, 12.0)
    try:
        # LE: ecliptic of date
        le_date = ephem.calc_ut(jd, 0, F | S)
        # LE: J2000
        le_j2k = ephem.calc_ut(jd, 0, F | S | J | N)

        prec_offset = (le_date[0][0] - le_j2k[0][0]) * 3600.0  # arcsec
        years_from_j2000 = year - 2000.0
        expected_offset = years_from_j2000 * 50.29  # ~50.29"/yr

        diff = abs(prec_offset - expected_offset)
        # Allow 5" per decade of nonlinear precession terms
        tol = 5.0 + abs(years_from_j2000) * 0.1
        if diff < tol:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL {year}: prec={prec_offset:.2f}" expected~{expected_offset:.2f}" diff={diff:.2f}"'
            )
    except Exception as e:
        errors += 1

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: SE vs LE precession offset agreement
# ============================================================
print("\n=== P2: SE vs LE precession offset ===")

for body, name in [(0, "Sun"), (4, "Mars"), (5, "Jupiter"), (6, "Saturn")]:
    for year in range(1960, 2041, 5):
        jd = swe.julday(year, 7, 1, 12.0)
        try:
            se_date = swe.calc_ut(jd, body, F | S)
            se_j2k = swe.calc_ut(jd, body, F | S | J | N)
            le_date = ephem.calc_ut(jd, body, F | S)
            le_j2k = ephem.calc_ut(jd, body, F | S | J | N)

            se_prec = (se_date[0][0] - se_j2k[0][0]) * 3600.0
            le_prec = (le_date[0][0] - le_j2k[0][0]) * 3600.0

            diff = abs(se_prec - le_prec)
            if diff < 1.0:  # Should agree within 1"
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(
                        f'  FAIL {name} {year}: SE_prec={se_prec:.2f}" LE_prec={le_prec:.2f}" diff={diff:.2f}"'
                    )
        except Exception as e:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Obliquity values SE vs LE
# ============================================================
print("\n=== P3: Obliquity (mean and true) ===")

for year in range(1800, 2201, 10):
    jd = swe.julday(year, 1, 1, 12.0)
    try:
        # Get obliquity via calc of ECL_NUT (body=-1 in some APIs)
        # Use nutation values from both
        se_nut = swe.calc_ut(jd, swe.ECL_NUT, 0)

        # LE uses internal nutation
        le_nut_result = ephem.calc_ut(jd, -1, 0)  # ECL_NUT = -1

        # Compare true obliquity (index 0)
        se_true_obl = se_nut[0][0]
        le_true_obl = le_nut_result[0][0]
        diff = abs(se_true_obl - le_true_obl) * 3600.0
        if diff < 1.0:
            passed += 1
        else:
            failed += 1
            if failed <= 5:
                print(
                    f'  FAIL obliquity {year}: SE={se_true_obl:.8f} LE={le_true_obl:.8f} diff={diff:.3f}"'
                )

        # Compare mean obliquity (index 1)
        se_mean_obl = se_nut[0][1]
        le_mean_obl = le_nut_result[0][1]
        diff2 = abs(se_mean_obl - le_mean_obl) * 3600.0
        if diff2 < 1.0:
            passed += 1
        else:
            failed += 1
    except Exception as e:
        errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Nutation in longitude and obliquity
# ============================================================
print("\n=== P4: Nutation components ===")

for year in range(1950, 2051, 2):
    jd = swe.julday(year, 1, 1, 12.0)
    try:
        se_nut = swe.calc_ut(jd, swe.ECL_NUT, 0)
        le_nut = ephem.calc_ut(jd, -1, 0)

        # Nutation in longitude (index 2)
        se_dpsi = se_nut[0][2]
        le_dpsi = le_nut[0][2]
        diff_dpsi = abs(se_dpsi - le_dpsi) * 3600.0
        if diff_dpsi < 0.5:
            passed += 1
        else:
            failed += 1
            if failed <= 5:
                print(
                    f'  FAIL dpsi {year}: SE={se_dpsi:.8f} LE={le_dpsi:.8f} diff={diff_dpsi:.4f}"'
                )

        # Nutation in obliquity (index 3)
        se_deps = se_nut[0][3]
        le_deps = le_nut[0][3]
        diff_deps = abs(se_deps - le_deps) * 3600.0
        if diff_deps < 0.5:
            passed += 1
        else:
            failed += 1
    except Exception as e:
        errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Sidereal time SE vs LE
# ============================================================
print("\n=== P5: Sidereal time ===")

for year in range(1950, 2051, 2):
    jd = swe.julday(year, 1, 1, 12.0)
    try:
        se_st = swe.sidtime(jd)
        le_st = ephem.sidtime(jd)
        diff_s = abs(se_st - le_st) * 3600.0  # hours to seconds
        if diff_s < 0.1:  # 0.1 second
            passed += 1
        else:
            failed += 1
            if failed <= 5:
                print(
                    f"  FAIL sidtime {year}: SE={se_st:.8f}h LE={le_st:.8f}h diff={diff_s:.4f}s"
                )
    except Exception as e:
        errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: Sidereal time at different UT hours
# ============================================================
print("\n=== P6: Sidereal time at different hours ===")

jd_base = swe.julday(2020, 3, 20, 0.0)  # Vernal equinox
for hour in range(0, 24, 3):
    jd = jd_base + hour / 24.0
    try:
        se_st = swe.sidtime(jd)
        le_st = ephem.sidtime(jd)
        diff_s = abs(se_st - le_st) * 3600.0
        if diff_s < 0.1:
            passed += 1
        else:
            failed += 1
    except Exception as e:
        errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 94 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
