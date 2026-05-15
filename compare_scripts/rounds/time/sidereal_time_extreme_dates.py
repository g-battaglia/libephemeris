#!/usr/bin/env python3
"""Round 52: Sidereal Time at Extreme Dates.

Tests sidtime() and sidtime0() across a wide range of dates,
comparing libephemeris against pyswisseph.

Phases:
  P1: sidtime at 50 dates spanning 1800-2200
  P2: sidtime0 with varying obliquity/nutation
  P3: Sidereal time consistency (sidtime vs sidtime0)
  P4: Deep historical dates (1000-1800 CE)
  P5: Sub-second precision at modern dates
  P6: Extreme future dates (2200-2500)
"""

from __future__ import annotations

import math
import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
results = {"passed": [], "failed": [], "errors": []}


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


def phase1():
    """sidtime at 50 dates spanning 1800-2200."""
    global errors
    print("\n=== P1: sidtime 1800-2200 ===")

    base_jd = 2378497.0  # 1800 Jan 1
    for i in range(50):
        jd = base_jd + i * 2922.0  # ~8 year steps
        year = 1800 + i * 8
        try:
            se_st = swe.sidtime(jd)
            le_st = ephem.sidtime(jd)

            # Sidereal time in hours (0-24)
            diff_h = abs(se_st - le_st)
            if diff_h > 12:
                diff_h = 24 - diff_h
            diff_s = diff_h * 3600  # to seconds

            # Tolerance: 0.5 seconds for modern era, relaxed for historical
            years_from_2000 = abs(year - 2000)
            tol = 0.5 + years_from_2000 * 0.005
            ok = diff_s < tol
            detail = f"SE={se_st:.8f}h LE={le_st:.8f}h diff={diff_s:.4f}s"
            record("P1", f"Y{year}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P1 Y{year}: {e}")


def phase2():
    """sidtime0 with varying obliquity/nutation."""
    global errors
    print("\n=== P2: sidtime0 ===")

    test_dates = [
        ("J2000", 2451545.0),
        ("2024Mar", 2460389.5),
        ("1950Jan", 2433283.0),
        ("1900Jan", 2415021.0),
        ("2100Jan", 2488070.0),
    ]

    # Test with different obliquity/nutation values
    obliquity_nutation_pairs = [
        (23.4393, 0.0, "eps23.44_nut0"),
        (23.4393, -0.00478, "eps23.44_nut-17"),
        (23.4393, 0.00478, "eps23.44_nut+17"),
        (23.0, 0.0, "eps23.0_nut0"),
        (24.0, 0.0, "eps24.0_nut0"),
    ]

    for date_name, jd in test_dates:
        for eps, nut, desc in obliquity_nutation_pairs:
            try:
                se_st0 = swe.sidtime0(jd, eps, nut)
                le_st0 = ephem.sidtime0(jd, eps, nut)

                diff_h = abs(se_st0 - le_st0)
                if diff_h > 12:
                    diff_h = 24 - diff_h
                diff_s = diff_h * 3600

                tol = 0.5
                ok = diff_s < tol
                detail = f"SE={se_st0:.8f}h LE={le_st0:.8f}h diff={diff_s:.4f}s"
                record("P2", f"{date_name} {desc}", ok, detail)

            except Exception as e:
                errors += 1
                results["errors"].append(f"P2 {date_name} {desc}: {e}")


def phase3():
    """Sidereal time consistency: sidtime should equal sidtime0 with correct eps/nut."""
    global errors
    print("\n=== P3: sidtime vs sidtime0 consistency ===")

    test_dates = [
        ("J2000", 2451545.0),
        ("2024Mar", 2460389.5),
        ("1980Jul", 2444407.0),
    ]

    for date_name, jd in test_dates:
        try:
            # Get sidereal time from both methods
            le_st = ephem.sidtime(jd)

            # Get obliquity and nutation for this date
            eps_true = ephem.calc_ut(jd, ephem.ECL_NUT, 0)[0]
            eps = eps_true[0]  # true obliquity
            nut = eps_true[2]  # nutation in longitude (degrees)

            le_st0 = ephem.sidtime0(jd, eps, nut)

            diff_h = abs(le_st - le_st0)
            if diff_h > 12:
                diff_h = 24 - diff_h
            diff_s = diff_h * 3600

            # Should be very close (both use same internal computation)
            tol = 0.01  # 10 ms
            ok = diff_s < tol
            detail = f"sidtime={le_st:.8f}h sidtime0={le_st0:.8f}h diff={diff_s:.6f}s"
            record("P3", f"consistency {date_name}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P3 {date_name}: {e}")


def phase4():
    """Deep historical dates (1000-1800 CE)."""
    global errors
    print("\n=== P4: Deep historical sidereal time (1000-1800) ===")

    for year in range(1000, 1800, 50):
        try:
            jd = ephem.julday(year, 1, 1, 12.0)

            se_st = swe.sidtime(jd)
            le_st = ephem.sidtime(jd)

            diff_h = abs(se_st - le_st)
            if diff_h > 12:
                diff_h = 24 - diff_h
            diff_s = diff_h * 3600

            # Historical: allow up to 2 seconds
            tol = 2.0
            ok = diff_s < tol
            detail = f"SE={se_st:.6f}h LE={le_st:.6f}h diff={diff_s:.3f}s"
            record("P4", f"Y{year}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 Y{year}: {e}")


def phase5():
    """Sub-second precision at modern dates."""
    global errors
    print("\n=== P5: Sub-second precision (modern) ===")

    # Test at hourly intervals for one day at J2000
    jd_base = 2451545.0
    for hour in range(24):
        jd = jd_base + hour / 24.0
        try:
            se_st = swe.sidtime(jd)
            le_st = ephem.sidtime(jd)

            diff_h = abs(se_st - le_st)
            if diff_h > 12:
                diff_h = 24 - diff_h
            diff_s = diff_h * 3600

            # Very tight: 0.1 seconds
            tol = 0.1
            ok = diff_s < tol
            detail = f"diff={diff_s:.6f}s"
            record("P5", f"J2000+{hour}h", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 {hour}h: {e}")

    # Also test at different dates with high precision
    precision_dates = [
        ("2024Jan01", 2460310.5),
        ("2024Jul01", 2460492.5),
        ("2010Jun15", 2455362.5),
        ("1999Dec31", 2451544.5),
    ]

    for name, jd in precision_dates:
        try:
            se_st = swe.sidtime(jd)
            le_st = ephem.sidtime(jd)

            diff_h = abs(se_st - le_st)
            if diff_h > 12:
                diff_h = 24 - diff_h
            diff_s = diff_h * 3600

            tol = 0.1
            ok = diff_s < tol
            detail = f"SE={se_st:.8f}h LE={le_st:.8f}h diff={diff_s:.6f}s"
            record("P5", f"precision {name}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 {name}: {e}")


def phase6():
    """Extreme future dates (2200-2500)."""
    global errors
    print("\n=== P6: Future sidereal time (2200-2500) ===")

    for year in range(2200, 2501, 25):
        try:
            jd = ephem.julday(year, 7, 1, 12.0)

            se_st = swe.sidtime(jd)
            le_st = ephem.sidtime(jd)

            diff_h = abs(se_st - le_st)
            if diff_h > 12:
                diff_h = 24 - diff_h
            diff_s = diff_h * 3600

            years_ahead = year - 2025
            tol = 1.0 + years_ahead * 0.005
            ok = diff_s < tol
            detail = f"SE={se_st:.6f}h LE={le_st:.6f}h diff={diff_s:.3f}s"
            record("P6", f"Y{year}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P6 Y{year}: {e}")


def main():
    print("=" * 70)
    print("ROUND 52: Sidereal Time at Extreme Dates")
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

    total = passed + failed + errors
    pct = 100 * passed / total if total else 0
    print("\n" + "=" * 70)
    print(f"ROUND 52 FINAL: {passed}/{total} passed ({pct:.1f}%)")
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
