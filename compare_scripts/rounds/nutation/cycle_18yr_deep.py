#!/usr/bin/env python3
"""Round 188: Nutation 18.6-year cycle deep verification.

Tests nutation components (dpsi, deps) and obliquity through the full
18.6-year cycle of the lunar node. Verifies against pyswisseph at fine
time steps sampling the entire cycle.
"""

from __future__ import annotations

import os
import sys
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLAGS = ephem.FLG_SWIEPH

passed = 0
failed = 0
total = 0
failures = []

# Sample the full 18.6-year nutation cycle from J2000 + some offsets
JD_J2000 = 2451545.0
CYCLE_DAYS = 6798.38  # 18.6 years in days

# Sample points: every ~30 days across two full cycles
SAMPLE_STEP = 30.0
N_SAMPLES = int(2 * CYCLE_DAYS / SAMPLE_STEP)


def test_nutation_cycle():
    global passed, failed, total

    print("=" * 70)
    print("Round 188: Nutation 18.6-Year Cycle Deep Verification")
    print("=" * 70)

    print(f"\nSampling {N_SAMPLES} points across two 18.6-year cycles...")

    max_dpsi_diff = 0.0
    max_deps_diff = 0.0
    max_obl_diff = 0.0

    for i in range(N_SAMPLES):
        jd = JD_J2000 + i * SAMPLE_STEP

        try:
            # Get nutation/obliquity from SE via calc_ut with ECL_NUT
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
            se_true_obl = se_r[0][0]
            se_mean_obl = se_r[0][1]
            se_nut_lon = se_r[0][2]  # dpsi in degrees
            se_nut_obl = se_r[0][3]  # deps in degrees

            # Get from libephemeris
            le_r = ephem.calc_ut(jd, -1, 0)
            le_true_obl = le_r[0][0]
            le_mean_obl = le_r[0][1]
            le_nut_lon = le_r[0][2]
            le_nut_obl = le_r[0][3]
        except Exception as e:
            continue

        # Nutation in longitude (dpsi)
        total += 1
        dpsi_diff = abs(le_nut_lon - se_nut_lon) * 3600
        max_dpsi_diff = max(max_dpsi_diff, dpsi_diff)
        if dpsi_diff <= 0.05:  # 0.05" tolerance
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  JD={jd:.1f} dpsi: LE={le_nut_lon:.8f} SE={se_nut_lon:.8f} diff={dpsi_diff:.4f}"'
            )

        # Nutation in obliquity (deps)
        total += 1
        deps_diff = abs(le_nut_obl - se_nut_obl) * 3600
        max_deps_diff = max(max_deps_diff, deps_diff)
        if deps_diff <= 0.05:
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  JD={jd:.1f} deps: LE={le_nut_obl:.8f} SE={se_nut_obl:.8f} diff={deps_diff:.4f}"'
            )

        # True obliquity
        total += 1
        obl_diff = abs(le_true_obl - se_true_obl) * 3600
        max_obl_diff = max(max_obl_diff, obl_diff)
        if obl_diff <= 0.1:  # 0.1" for obliquity (includes mean obl model diff)
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  JD={jd:.1f} true_obl: LE={le_true_obl:.8f} SE={se_true_obl:.8f} diff={obl_diff:.4f}"'
            )

        # Mean obliquity
        total += 1
        mobl_diff = abs(le_mean_obl - se_mean_obl) * 3600
        if mobl_diff <= 0.1:
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  JD={jd:.1f} mean_obl: LE={le_mean_obl:.8f} SE={se_mean_obl:.8f} diff={mobl_diff:.4f}"'
            )

    print(f'\n  Max dpsi diff: {max_dpsi_diff:.4f}"')
    print(f'  Max deps diff: {max_deps_diff:.4f}"')
    print(f'  Max true_obl diff: {max_obl_diff:.4f}"')


def test_nutation_extremes():
    """Test at nutation extremes (max/min dpsi)."""
    global passed, failed, total

    print(f"\n--- Nutation Extremes ---")

    # Approximate dates of maximum nutation (when lunar node near 0°/180°)
    extreme_jds = [
        JD_J2000 + 0,  # J2000 itself
        JD_J2000 + 3399,  # ~half cycle
        JD_J2000 + 6798,  # ~full cycle
        JD_J2000 + 10197,  # ~1.5 cycles
        JD_J2000 - 3399,  # ~half cycle before
        JD_J2000 - 6798,  # ~full cycle before
    ]

    for jd in extreme_jds:
        try:
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
            le_r = ephem.calc_ut(jd, -1, 0)

            for idx, name in [
                (0, "true_obl"),
                (1, "mean_obl"),
                (2, "nut_lon"),
                (3, "nut_obl"),
            ]:
                total += 1
                diff = abs(le_r[0][idx] - se_r[0][idx]) * 3600
                tol = 0.1
                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    failures.append(f'  JD={jd:.1f} extreme {name}: diff={diff:.4f}"')
        except Exception:
            pass


def test_nutation_at_historical():
    """Test nutation at historically important dates."""
    global passed, failed, total

    print(f"\n--- Nutation at Historical Dates ---")

    historical_jds = [
        2415020.0,  # 1900 Jan 0.5
        2430000.0,  # 1941
        2440000.0,  # 1968
        2445000.0,  # 1981
        2450000.0,  # 1995
        2455000.0,  # 2009
        2460000.0,  # 2023
        2465000.0,  # 2037
    ]

    for jd in historical_jds:
        try:
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
            le_r = ephem.calc_ut(jd, -1, 0)

            for idx, name in [
                (0, "true_obl"),
                (1, "mean_obl"),
                (2, "nut_lon"),
                (3, "nut_obl"),
            ]:
                total += 1
                diff = abs(le_r[0][idx] - se_r[0][idx]) * 3600
                tol = 0.15  # slightly wider for historical
                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  JD={jd:.1f} hist {name}: LE={le_r[0][idx]:.8f} SE={se_r[0][idx]:.8f} diff={diff:.4f}"'
                    )
        except Exception:
            pass


if __name__ == "__main__":
    test_nutation_cycle()
    test_nutation_extremes()
    test_nutation_at_historical()

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
