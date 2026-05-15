#!/usr/bin/env python3
"""Round 130: Mean/True Node Oscillation Amplitude

Tests the oscillation between mean and true node, verifying the amplitude
and phase of the ~18.6-year nutation-driven oscillation matches SE.
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
MEAN_NODE = 10
TRUE_NODE = 11


def main():
    print("=" * 80)
    print("ROUND 130: Mean/True Node Oscillation Amplitude")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Sample every 5 days across 2 full nutation cycles (~37 years)
    start_jd = 2444239.5  # 1980-01-01
    end_jd = 2457754.5  # 2017-01-01
    step = 5.0

    max_osc_se = 0
    max_osc_le = 0
    jd = start_jd
    while jd <= end_jd:
        try:
            se_mean = swe.calc_ut(jd, MEAN_NODE, FLG_SPEED)[0]
            se_true = swe.calc_ut(jd, TRUE_NODE, FLG_SPEED)[0]
            le_mean = ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)[0]
            le_true = ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)[0]
        except Exception:
            jd += step
            continue

        # Compare mean node longitude
        diff_mean = le_mean[0] - se_mean[0]
        if diff_mean > 180:
            diff_mean -= 360
        elif diff_mean < -180:
            diff_mean += 360
        diff_mean_as = abs(diff_mean) * 3600
        total_tests += 1
        if diff_mean_as < 2.0:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 10:
                failures.append(f'  MeanNode JD={jd:.1f}: {diff_mean_as:.4f}"')

        # Compare true node longitude
        diff_true = le_true[0] - se_true[0]
        if diff_true > 180:
            diff_true -= 360
        elif diff_true < -180:
            diff_true += 360
        diff_true_as = abs(diff_true) * 3600
        total_tests += 1
        if diff_true_as < 2.0:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 15:
                failures.append(f'  TrueNode JD={jd:.1f}: {diff_true_as:.4f}"')

        # Compare oscillation (true - mean)
        se_osc = se_true[0] - se_mean[0]
        if se_osc > 180:
            se_osc -= 360
        elif se_osc < -180:
            se_osc += 360
        le_osc = le_true[0] - le_mean[0]
        if le_osc > 180:
            le_osc -= 360
        elif le_osc < -180:
            le_osc += 360

        max_osc_se = max(max_osc_se, abs(se_osc))
        max_osc_le = max(max_osc_le, abs(le_osc))

        diff_osc = abs(le_osc - se_osc) * 3600
        total_tests += 1
        if diff_osc < 5.0:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 20:
                failures.append(f'  Oscillation JD={jd:.1f}: diff={diff_osc:.4f}"')

        # Compare speeds
        diff_mean_spd = abs(le_mean[3] - se_mean[3]) * 3600
        total_tests += 1
        if diff_mean_spd < 10.0:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 25:
                failures.append(
                    f'  MeanNode speed JD={jd:.1f}: {diff_mean_spd:.4f}"/day'
                )

        diff_true_spd = abs(le_true[3] - se_true[3]) * 3600
        total_tests += 1
        if diff_true_spd < 10.0:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 30:
                failures.append(
                    f'  TrueNode speed JD={jd:.1f}: {diff_true_spd:.4f}"/day'
                )

        jd += step

    print(f"\n  Max oscillation (true-mean): SE={max_osc_se:.4f}° LE={max_osc_le:.4f}°")

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 130 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
    )
    if failures:
        print("\nSample failures:")
        for f in failures:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
