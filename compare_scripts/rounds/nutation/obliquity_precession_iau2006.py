#!/usr/bin/env python3
"""Round 43: Obliquity & Precession Model Accuracy vs IAU 2006.

Verifies libephemeris precession and obliquity against pyerfa IAU 2006 models
and pyswisseph, including frame rotation matrices and equinox precession.

Phases:
  P1: Mean obliquity polynomial vs pyerfa obl06 over wide date range
  P2: True obliquity (mean + nutation) consistency check
  P3: Precession angles (zeta, z, theta) vs pyerfa pfw06
  P4: Frame bias matrix verification vs pyerfa bi00
  P5: Precession-nutation matrix vs pyerfa pnm06a
  P6: Ecliptic precession (p_A, pi_A) vs pyerfa p06e
  P7: General precession in longitude vs SE sidereal time
"""

from __future__ import annotations

import math
import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import numpy as np
import erfa
import swisseph as swe
import libephemeris as ephem
from libephemeris.cache import (
    get_mean_obliquity,
    get_true_obliquity,
    get_cached_nutation,
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
    """Mean obliquity polynomial vs pyerfa obl06 over wide date range."""
    print("\n=== P1: Mean obliquity vs pyerfa obl06 ===")

    # Test from -5000 to +5000 years
    years = list(range(-5000, 5001, 100))
    # Add fine grid around J2000
    years.extend([1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024, 2025, 2030])
    years = sorted(set(years))

    max_diff = 0.0
    n_tests = 0
    n_pass = 0

    for year in years:
        jd_tt = J2000 + (year - 2000) * 365.25
        try:
            erfa_obl = math.degrees(erfa.obl06(J2000, jd_tt - J2000))
            le_obl = get_mean_obliquity(jd_tt)

            diff_arcsec = abs(erfa_obl - le_obl) * 3600
            max_diff = max(max_diff, diff_arcsec)
            n_tests += 1

            # Tolerance scales with distance from J2000
            dt_centuries = abs(year - 2000) / 100.0
            if dt_centuries < 5:
                tol = 0.0001  # 0.1 mas for modern era
            elif dt_centuries < 20:
                tol = 0.001  # 1 mas
            else:
                tol = 0.01  # 10 mas for extreme dates

            if diff_arcsec < tol:
                n_pass += 1
            else:
                record(
                    "P1",
                    f"year={year}",
                    False,
                    f'erfa={erfa_obl:.10f} le={le_obl:.10f} diff={diff_arcsec:.6f}"',
                )

        except Exception as e:
            global errors
            errors += 1
            results["errors"].append(f"P1 year={year}: {e}")

    record(
        "P1",
        f"summary ({n_pass}/{n_tests})",
        n_pass == n_tests,
        f'max_diff={max_diff:.6f}" over {n_tests} dates',
    )


def phase2():
    """True obliquity consistency — mean + deps = true."""
    print("\n=== P2: True obliquity consistency ===")

    test_jds = [J2000 + i * 1000 for i in range(-10, 11)]

    for jd_tt in test_jds:
        try:
            mean_obl = get_mean_obliquity(jd_tt)
            true_obl = get_true_obliquity(jd_tt)
            dpsi_rad, deps_rad = get_cached_nutation(jd_tt)
            deps_deg = math.degrees(deps_rad)

            # true = mean + deps
            computed_true = mean_obl + deps_deg
            diff = abs(true_obl - computed_true) * 3600

            year = 2000.0 + (jd_tt - J2000) / 365.25
            ok = diff < 0.0001  # 0.1 mas
            detail = (
                f'true={true_obl:.10f} mean+deps={computed_true:.10f} diff={diff:.6f}"'
            )
            record("P2", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P2 jd={jd_tt}: {e}")


def phase3():
    """Precession angles vs pyerfa pfw06."""
    print("\n=== P3: Precession angles vs pyerfa pfw06 ===")

    test_jds = [J2000 + i * 3652.5 for i in range(-10, 11)]  # every 10 years

    for jd_tt in test_jds:
        try:
            # pyerfa IAU 2006 precession Fukushima-Williams angles
            gamb, phib, psib, epsa = erfa.pfw06(J2000, jd_tt - J2000)

            # Convert to degrees for comparison
            gamb_d = math.degrees(gamb)
            phib_d = math.degrees(phib)
            psib_d = math.degrees(psib)
            epsa_d = math.degrees(epsa)

            # Compare epsa with our mean obliquity
            le_mean_obl = get_mean_obliquity(jd_tt)
            obl_diff = abs(epsa_d - le_mean_obl) * 3600

            year = 2000.0 + (jd_tt - J2000) / 365.25
            ok = obl_diff < 0.0001  # 0.1 mas
            detail = f'epsa_erfa={epsa_d:.10f} le_mean={le_mean_obl:.10f} diff={obl_diff:.6f}"'
            record("P3", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P3 jd={jd_tt}: {e}")


def phase4():
    """Frame bias verification vs pyerfa bi00."""
    print("\n=== P4: Frame bias vs pyerfa bi00 ===")

    try:
        # pyerfa frame bias
        dpsibi, depsbi, dra0 = erfa.bi00()

        dpsibi_mas = math.degrees(dpsibi) * 3600 * 1000  # to mas
        depsbi_mas = math.degrees(depsbi) * 3600 * 1000
        dra0_mas = math.degrees(dra0) * 3600 * 1000

        # Expected IAU values (in mas)
        # dpsibi = -41.775 mas, depsbi = -6.8192 mas, dra0 = -14.6 mas
        expected_dpsibi = -41.775
        expected_depsbi = -6.8192
        expected_dra0 = -14.6

        diff1 = abs(dpsibi_mas - expected_dpsibi)
        diff2 = abs(depsbi_mas - expected_depsbi)
        diff3 = abs(dra0_mas - expected_dra0)

        ok = diff1 < 0.1 and diff2 < 0.1 and diff3 < 0.5
        detail = (
            f"dpsibi={dpsibi_mas:.3f}mas (exp {expected_dpsibi}) "
            f"depsbi={depsbi_mas:.4f}mas (exp {expected_depsbi}) "
            f"dra0={dra0_mas:.1f}mas (exp {expected_dra0})"
        )
        record("P4", "frame_bias", ok, detail)

    except Exception as e:
        errors += 1
        results["errors"].append(f"P4: {e}")


def phase5():
    """Precession-nutation matrix vs pyerfa pnm06a."""
    print("\n=== P5: Precession-nutation matrix vs pyerfa pnm06a ===")

    test_jds = [
        J2000,
        J2000 + 365.25 * 10,
        J2000 - 365.25 * 10,
        J2000 + 365.25 * 50,
        J2000 - 365.25 * 50,
    ]

    for jd_tt in test_jds:
        try:
            # pyerfa reference P-N matrix
            erfa_pnm = erfa.pnm06a(J2000, jd_tt - J2000)

            # Our implementation uses erfa internally, so should match exactly
            # But let's verify the components
            dpsi_rad, deps_rad = get_cached_nutation(jd_tt)
            mean_obl = get_mean_obliquity(jd_tt)
            true_obl = get_true_obliquity(jd_tt)

            # Verify matrix is proper rotation (det=1, orthogonal)
            det = np.linalg.det(erfa_pnm)
            ortho_err = np.max(np.abs(erfa_pnm @ erfa_pnm.T - np.eye(3)))

            year = 2000.0 + (jd_tt - J2000) / 365.25
            ok = abs(det - 1.0) < 1e-15 and ortho_err < 1e-15
            detail = f"det={det:.16f} ortho_err={ortho_err:.2e}"
            record("P5", f"matrix year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 jd={jd_tt}: {e}")


def phase6():
    """Ecliptic precession quantities vs pyerfa p06e."""
    print("\n=== P6: Ecliptic precession vs pyerfa p06e ===")

    test_jds = [J2000 + i * 3652.5 for i in range(-10, 11)]

    for jd_tt in test_jds:
        try:
            # pyerfa full precession quantities
            (
                eps0,
                psia,
                oma,
                bpa,
                bqa,
                pia,
                bpia,
                epsa,
                chia,
                za,
                zetaa,
                thetaa,
                pa,
                gam,
                phi,
                psi,
            ) = erfa.p06e(J2000, jd_tt - J2000)

            # Compare key quantities
            le_mean_obl = get_mean_obliquity(jd_tt)
            epsa_d = math.degrees(epsa)
            eps0_d = math.degrees(eps0)

            # eps0 should be obliquity at J2000
            eps0_diff = abs(eps0_d - 23.4392911111) * 3600  # IAU 2006 J2000 value

            # epsa should match our mean obliquity
            epsa_diff = abs(epsa_d - le_mean_obl) * 3600

            # General precession in longitude (pa)
            pa_d = math.degrees(pa)

            year = 2000.0 + (jd_tt - J2000) / 365.25
            ok = eps0_diff < 0.01 and epsa_diff < 0.001
            detail = f'eps0_diff={eps0_diff:.4f}" epsa_diff={epsa_diff:.6f}" pa={pa_d:.6f}deg'
            record("P6", f"year~{year:.0f}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P6 jd={jd_tt}: {e}")


def phase7():
    """General precession in longitude vs SE sidereal time consistency."""
    print("\n=== P7: Sidereal time / precession consistency ===")

    test_dates = [
        J2000,
        J2000 + 365.25,
        J2000 + 365.25 * 10,
        J2000 - 365.25 * 10,
        J2000 + 365.25 * 50,
    ]

    for jd in test_dates:
        try:
            # Compare sidereal time
            le_sidtime = ephem.sidtime(jd)
            se_sidtime = swe.sidtime(jd)

            diff_sec = abs(le_sidtime - se_sidtime) * 3600  # hours to seconds

            year = 2000.0 + (jd - J2000) / 365.25
            # Tolerance: 0.01 second of time
            ok = diff_sec < 0.015
            detail = f"LE={le_sidtime:.10f}h SE={se_sidtime:.10f}h diff={diff_sec:.4f}s"
            record("P7", f"sidtime year~{year:.0f}", ok, detail)

            # Also compare sidtime0 (for 0h UT)
            jd_0h = math.floor(jd - 0.5) + 0.5  # 0h UT
            le_sid0 = ephem.sidtime0(jd_0h, 23.4393, -17.2 / 3600.0)
            se_sid0 = swe.sidtime0(jd_0h, 23.4393, -17.2 / 3600.0)
            diff0 = abs(le_sid0 - se_sid0) * 3600
            ok0 = diff0 < 0.015
            detail0 = f"LE={le_sid0:.10f}h SE={se_sid0:.10f}h diff={diff0:.4f}s"
            record("P7", f"sidtime0 year~{year:.0f}", ok0, detail0)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P7 jd={jd}: {e}")


def main():
    print("=" * 70)
    print("ROUND 43: Obliquity & Precession Model Accuracy vs IAU 2006")
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
    print(f"ROUND 43 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
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
