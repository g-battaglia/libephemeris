#!/usr/bin/env python3
"""Round 124: Nutation Components at 18.6yr Cycle Nodes

Tests nutation (dpsi, deps) at critical points in the 18.6-year lunar nodal cycle
where nutation amplitude is at maximum/minimum.
"""

from __future__ import annotations
import os, sys, math

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256


def main():
    print("=" * 80)
    print("ROUND 124: Nutation Components at 18.6yr Cycle Nodes")
    print("=" * 80)

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Sample JDs spanning ~4 full nutation cycles (1940-2050)
    # Step every 30 days for fine coverage
    start_jd = 2429630.0  # ~1940
    end_jd = 2469808.0  # ~2050
    step = 30.0

    jd = start_jd
    while jd <= end_jd:
        # Get nutation from both
        try:
            se_nut = swe.calc_ut(jd, -1, 0)  # ECL_NUT = -1
            # Returns (true_obl, mean_obl, dpsi, deps, 0, 0)
            se_true_obl = se_nut[0][0]
            se_mean_obl = se_nut[0][1]
            se_dpsi = se_nut[0][2]
            se_deps = se_nut[0][3]
        except Exception:
            jd += step
            continue

        try:
            le_nut = ephem.calc_ut(jd, -1, 0)
            le_true_obl = le_nut[0][0]
            le_mean_obl = le_nut[0][1]
            le_dpsi = le_nut[0][2]
            le_deps = le_nut[0][3]
        except Exception:
            jd += step
            continue

        # Compare true obliquity
        diff_true = abs(le_true_obl - se_true_obl) * 3600
        total_tests += 1
        if diff_true < 0.5:  # 0.5" for obliquity
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 10:
                failures.append(f'  JD={jd:.1f}: true_obl diff={diff_true:.4f}"')

        # Compare mean obliquity
        diff_mean = abs(le_mean_obl - se_mean_obl) * 3600
        total_tests += 1
        if diff_mean < 0.5:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 15:
                failures.append(f'  JD={jd:.1f}: mean_obl diff={diff_mean:.4f}"')

        # Compare nutation in longitude (dpsi)
        diff_dpsi = abs(le_dpsi - se_dpsi) * 3600
        total_tests += 1
        if diff_dpsi < 0.5:  # 0.5" tolerance
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 20:
                failures.append(f'  JD={jd:.1f}: dpsi diff={diff_dpsi:.4f}"')

        # Compare nutation in obliquity (deps)
        diff_deps = abs(le_deps - se_deps) * 3600
        total_tests += 1
        if diff_deps < 0.5:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 25:
                failures.append(f'  JD={jd:.1f}: deps diff={diff_deps:.4f}"')

        jd += step

    print(f"\n{'=' * 80}")
    pct = 100 * total_pass / total_tests if total_tests else 0
    print(f"ROUND 124 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
    print(f"  Failures: {total_fail}")
    print("=" * 80)
    if failures:
        print("\nSample failures:")
        for f in failures:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
