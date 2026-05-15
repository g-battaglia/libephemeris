#!/usr/bin/env python3
"""Round 42: Nutation Components IAU 2000B Deep Verification.

Verifies libephemeris nutation/obliquity against pyswisseph and pyerfa independently.

Phases:
  P1: ECL_NUT query — nutation angles and obliquity vs pyswisseph
  P2: Nutation angles vs pyerfa IAU 2006/2000A directly
  P3: Mean obliquity vs pyerfa IAU 2006 directly
  P4: Nutation at historical/future dates (1000 BCE to 3000 CE)
  P5: Nutation consistency — dpsi/deps match between calc_ut and cache
  P6: FLG_NONUT flag — verify mean vs true obliquity selection
  P7: Nutation speed — dpsi/deps rate of change verification
"""

from __future__ import annotations

import math
import os
import sys
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import erfa
import swisseph as swe
import libephemeris as ephem
from libephemeris.cache import (
    get_cached_nutation,
    get_nutation_degrees,
    get_mean_obliquity,
    get_true_obliquity,
)

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
results = {"passed": [], "failed": [], "errors": []}

J2000 = 2451545.0


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


def phase1():
    """ECL_NUT query — compare LE vs SE for nutation and obliquity."""
    print("\n=== P1: ECL_NUT nutation/obliquity vs pyswisseph ===")

    test_dates = {
        "J2000": 2451545.0,
        "2024.0": 2460310.5,
        "2024.5": 2460493.0,
        "1990.0": 2447892.5,
        "1950.0": 2433282.5,
        "1900.0": 2415020.5,
        "2050.0": 2469807.5,
        "2100.0": 2488069.5,
        "1800.0": 2378496.5,
        "1600.0": 2305447.5,
    }

    for label, jd in test_dates.items():
        try:
            se_result = swe.calc_ut(jd, swe.ECL_NUT)[0]
            le_result = ephem.calc_ut(jd, ephem.ECL_NUT, 0)[0]

            # SE returns: (true_obl, mean_obl, dpsi, deps, 0, 0)
            se_true_obl = se_result[0]
            se_mean_obl = se_result[1]
            se_dpsi = se_result[2]
            se_deps = se_result[3]

            le_true_obl = le_result[0]
            le_mean_obl = le_result[1]
            le_dpsi = le_result[2]
            le_deps = le_result[3]

            # Tolerances
            obl_tol = 0.1 / 3600.0  # 0.1 arcsecond
            nut_tol = 0.01 / 3600.0  # 0.01 arcsecond (10 mas)

            true_obl_diff = abs(se_true_obl - le_true_obl) * 3600
            mean_obl_diff = abs(se_mean_obl - le_mean_obl) * 3600
            dpsi_diff = abs(se_dpsi - le_dpsi) * 3600
            deps_diff = abs(se_deps - le_deps) * 3600

            all_ok = (
                true_obl_diff < obl_tol * 3600
                and mean_obl_diff < obl_tol * 3600
                and dpsi_diff < nut_tol * 3600
                and deps_diff < nut_tol * 3600
            )

            detail = (
                f'true_obl={true_obl_diff:.4f}" mean_obl={mean_obl_diff:.4f}" '
                f'dpsi={dpsi_diff:.4f}" deps={deps_diff:.4f}"'
            )
            record("P1", label, all_ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P1 {label}: {e}")


def phase2():
    """Nutation angles vs pyerfa IAU 2006/2000A directly."""
    print("\n=== P2: Nutation angles vs pyerfa nut06a ===")

    test_jds_tt = [
        J2000,
        J2000 + 365.25 * 24,  # 2024
        J2000 - 365.25 * 50,  # 1950
        J2000 + 365.25 * 50,  # 2050
        J2000 + 365.25 * 100,  # 2100
        J2000 - 365.25 * 100,  # 1900
        J2000 + 365.25 * 200,  # 2200
        J2000 - 365.25 * 200,  # 1800
    ]

    for jd_tt in test_jds_tt:
        try:
            # pyerfa reference
            erfa_dpsi, erfa_deps = erfa.nut06a(J2000, jd_tt - J2000)

            # libephemeris cached nutation (returns radians)
            le_dpsi_rad, le_deps_rad = get_cached_nutation(jd_tt)

            dpsi_diff_mas = abs(erfa_dpsi - le_dpsi_rad) * 206264806.0  # rad to mas
            deps_diff_mas = abs(erfa_deps - le_deps_rad) * 206264806.0

            year = 2000.0 + (jd_tt - J2000) / 365.25
            tol_mas = 0.001  # 1 microarcsecond tolerance (should be exact match)

            ok = dpsi_diff_mas < tol_mas and deps_diff_mas < tol_mas
            detail = (
                f"dpsi_diff={dpsi_diff_mas:.6f}mas deps_diff={deps_diff_mas:.6f}mas"
            )
            record("P2", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P2 jd={jd_tt}: {e}")


def phase3():
    """Mean obliquity vs pyerfa IAU 2006 directly."""
    print("\n=== P3: Mean obliquity vs pyerfa obl06 ===")

    test_jds_tt = [
        J2000,
        J2000 + 365.25 * 24,
        J2000 - 365.25 * 50,
        J2000 + 365.25 * 50,
        J2000 + 365.25 * 100,
        J2000 - 365.25 * 100,
        J2000 + 365.25 * 500,
        J2000 - 365.25 * 500,
    ]

    for jd_tt in test_jds_tt:
        try:
            erfa_obl = math.degrees(erfa.obl06(J2000, jd_tt - J2000))
            le_obl = get_mean_obliquity(jd_tt)

            diff_arcsec = abs(erfa_obl - le_obl) * 3600
            year = 2000.0 + (jd_tt - J2000) / 365.25
            tol = 0.0001  # 0.1 mas

            ok = diff_arcsec < tol
            detail = f'erfa={erfa_obl:.10f} le={le_obl:.10f} diff={diff_arcsec:.6f}"'
            record("P3", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P3 jd={jd_tt}: {e}")


def phase4():
    """Nutation at historical/future dates — wide date sweep."""
    print("\n=== P4: Nutation wide-date sweep vs SE ===")

    # Sweep from 1000 BCE to 3000 CE
    years = [
        -1000,
        -500,
        0,
        500,
        1000,
        1200,
        1400,
        1500,
        1600,
        1700,
        1800,
        1850,
        1900,
        1920,
        1940,
        1960,
        1980,
        1990,
        2000,
        2010,
        2020,
        2024,
        2025,
        2030,
        2050,
        2100,
        2200,
        2300,
        2500,
        3000,
    ]

    for year in years:
        jd = J2000 + (year - 2000) * 365.25
        try:
            se_result = swe.calc_ut(jd, swe.ECL_NUT)[0]
            le_result = ephem.calc_ut(jd, ephem.ECL_NUT, 0)[0]

            true_obl_diff = abs(se_result[0] - le_result[0]) * 3600
            mean_obl_diff = abs(se_result[1] - le_result[1]) * 3600
            dpsi_diff = abs(se_result[2] - le_result[2]) * 3600
            deps_diff = abs(se_result[3] - le_result[3]) * 3600

            # Tolerance: 0.5" for historical, 0.01" for modern
            if abs(year - 2000) < 200:
                tol = 0.01  # 10 mas
            elif abs(year - 2000) < 500:
                tol = 0.05  # 50 mas
            else:
                tol = 0.5  # 500 mas for ancient/far future

            all_ok = (
                true_obl_diff < tol
                and mean_obl_diff < tol
                and dpsi_diff < tol
                and deps_diff < tol
            )

            detail = (
                f'true_obl={true_obl_diff:.4f}" mean_obl={mean_obl_diff:.4f}" '
                f'dpsi={dpsi_diff:.4f}" deps={deps_diff:.4f}"'
            )
            record("P4", f"year={year}", all_ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 year={year}: {e}")


def phase5():
    """Internal consistency — nutation from calc_ut matches cache."""
    print("\n=== P5: Internal nutation consistency ===")

    test_jds = [
        J2000,
        J2000 + 100,
        J2000 + 1000,
        J2000 - 1000,
        J2000 + 10000,
        J2000 - 10000,
    ]

    for jd in test_jds:
        try:
            # Get from calc_ut
            le_result = ephem.calc_ut(jd, ephem.ECL_NUT, 0)[0]
            le_dpsi_deg = le_result[2]
            le_deps_deg = le_result[3]
            le_mean_obl = le_result[1]
            le_true_obl = le_result[0]

            # Get from cache
            from libephemeris.state import get_timescale

            ts = get_timescale()
            t = ts.ut1_jd(jd)
            jd_tt = t.tt

            cache_dpsi_deg, cache_deps_deg = get_nutation_degrees(jd_tt)
            cache_mean_obl = get_mean_obliquity(jd_tt)
            cache_true_obl = get_true_obliquity(jd_tt)

            dpsi_diff = abs(le_dpsi_deg - cache_dpsi_deg) * 3600
            deps_diff = abs(le_deps_deg - cache_deps_deg) * 3600
            mean_diff = abs(le_mean_obl - cache_mean_obl) * 3600
            true_diff = abs(le_true_obl - cache_true_obl) * 3600

            tol = 0.0001  # 0.1 mas
            ok = (
                dpsi_diff < tol
                and deps_diff < tol
                and mean_diff < tol
                and true_diff < tol
            )

            year = 2000.0 + (jd - J2000) / 365.25
            detail = (
                f'dpsi={dpsi_diff:.6f}" deps={deps_diff:.6f}" '
                f'mean={mean_diff:.6f}" true={true_diff:.6f}"'
            )
            record("P5", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 jd={jd}: {e}")


def phase6():
    """FLG_NONUT flag — mean vs true obliquity."""
    print("\n=== P6: FLG_NONUT flag verification ===")

    test_dates = [J2000, J2000 + 5000, J2000 - 5000]

    for jd in test_dates:
        try:
            # Without NONUT — should use true equinox
            le_pos, _ = ephem.calc_ut(jd, 0, ephem.FLG_SPEED)
            # With NONUT — should use mean equinox
            le_pos_nn, _ = ephem.calc_ut(
                jd, 0, ephem.FLG_SPEED | ephem.FLG_NONUT
            )

            # The positions should differ by approximately dpsi * cos(eps)
            le_nut = ephem.calc_ut(jd, ephem.ECL_NUT, 0)[0]
            dpsi = le_nut[2]  # degrees
            eps = le_nut[0]  # true obliquity

            lon_diff = le_pos[0] - le_pos_nn[0]
            expected_diff = dpsi  # nutation in longitude

            # The difference should be close to dpsi
            residual = abs(lon_diff - expected_diff) * 3600
            tol = 1.0  # 1 arcsecond tolerance (not exact due to other effects)

            year = 2000.0 + (jd - J2000) / 365.25
            ok = residual < tol
            detail = (
                f'lon_diff={lon_diff * 3600:.4f}" dpsi={dpsi * 3600:.4f}" '
                f'residual={residual:.4f}"'
            )
            record("P6", f"Sun year~{year:.0f}", ok, detail)

            # Also check SE does the same
            se_pos = swe.calc_ut(jd, 0, swe.FLG_SPEED)[0]
            se_pos_nn = swe.calc_ut(jd, 0, swe.FLG_SPEED | swe.FLG_NONUT)[0]
            se_lon_diff = se_pos[0] - se_pos_nn[0]

            cross_diff = abs(lon_diff - se_lon_diff) * 3600
            ok2 = cross_diff < 0.01  # 10 mas
            detail2 = f'LE_diff={lon_diff * 3600:.4f}" SE_diff={se_lon_diff * 3600:.4f}" cross={cross_diff:.4f}"'
            record("P6", f"Sun cross year~{year:.0f}", ok2, detail2)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P6 jd={jd}: {e}")


def phase7():
    """Nutation rate of change — finite-difference dpsi/deps speed."""
    print("\n=== P7: Nutation rate of change (speed) ===")

    test_jds = [J2000, J2000 + 5000, J2000 - 5000, J2000 + 15000]
    dt = 0.001  # 0.001 day ~ 86 seconds

    for jd in test_jds:
        try:
            # LE nutation at jd-dt, jd, jd+dt
            le_m = ephem.calc_ut(jd - dt, ephem.ECL_NUT, 0)[0]
            le_0 = ephem.calc_ut(jd, ephem.ECL_NUT, 0)[0]
            le_p = ephem.calc_ut(jd + dt, ephem.ECL_NUT, 0)[0]

            # SE nutation at jd-dt, jd, jd+dt
            se_m = swe.calc_ut(jd - dt, swe.ECL_NUT)[0]
            se_0 = swe.calc_ut(jd, swe.ECL_NUT)[0]
            se_p = swe.calc_ut(jd + dt, swe.ECL_NUT)[0]

            # Finite-difference dpsi speed
            le_dpsi_speed = (le_p[2] - le_m[2]) / (2 * dt)  # deg/day
            se_dpsi_speed = (se_p[2] - se_m[2]) / (2 * dt)

            le_deps_speed = (le_p[3] - le_m[3]) / (2 * dt)
            se_deps_speed = (se_p[3] - se_m[3]) / (2 * dt)

            dpsi_speed_diff = abs(le_dpsi_speed - se_dpsi_speed) * 3600  # arcsec/day
            deps_speed_diff = abs(le_deps_speed - se_deps_speed) * 3600

            tol = 0.001  # 1 mas/day
            year = 2000.0 + (jd - J2000) / 365.25
            ok = dpsi_speed_diff < tol and deps_speed_diff < tol
            detail = (
                f'dpsi_speed_diff={dpsi_speed_diff:.6f}"/day '
                f'deps_speed_diff={deps_speed_diff:.6f}"/day'
            )
            record("P7", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P7 jd={jd}: {e}")


def main():
    print("=" * 70)
    print("ROUND 42: Nutation Components IAU 2000B Deep Verification")
    print("=" * 70)

    phase1()
    print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

    phase2()
    print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

    phase3()
    print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

    phase4()
    print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

    phase5()
    print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

    phase6()
    print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

    phase7()
    print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

    total = passed + failed + errors
    print("\n" + "=" * 70)
    print(f"ROUND 42 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print("=" * 70)

    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:30]:
            print(f)

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:10]:
            print(e)

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
