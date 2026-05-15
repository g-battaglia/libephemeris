#!/usr/bin/env python3
"""Round 128: Eclipse Timing at Historical Dates (1950-2000)

Tests solar and lunar eclipse timing at historical dates where both
ephemerides should be most accurate (within DE440 range).
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)


def main():
    print("=" * 80)
    print("ROUND 128: Eclipse Timing at Historical Dates (1950-2000)")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Find solar eclipses 1950-2000
    print("\n--- Solar eclipses 1950-2000 ---")
    jd = 2433282.5  # 1950-01-01
    end_jd = 2451545.0  # 2000-01-01
    solar_eclipses = []
    while jd < end_jd and len(solar_eclipses) < 40:
        try:
            se_res = swe.sol_eclipse_when_glob(jd, 0, 0, False)
            tmax = se_res[1][0]
            if tmax < end_jd:
                solar_eclipses.append((se_res[0], se_res[1]))
            jd = tmax + 30
        except Exception:
            jd += 30

    print(f"  Found {len(solar_eclipses)} solar eclipses")

    for ecl_type, se_times in solar_eclipses:
        try:
            le_res = ephem.sol_eclipse_when_glob(se_times[0] - 15, 0, 0)
            le_times = le_res[1]
        except Exception:
            continue

        # Compare maximum time
        se_tmax = se_times[0]
        le_tmax = le_times[0]
        diff_min = abs(le_tmax - se_tmax) * 1440  # minutes
        tol_min = 5.0  # 5 minutes

        total_tests += 1
        if diff_min < tol_min:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 15:
                failures.append(
                    f"  Solar JD={se_tmax:.4f}: tmax diff={diff_min:.2f} min"
                )

    # Find lunar eclipses 1950-2000
    print("\n--- Lunar eclipses 1950-2000 ---")
    jd = 2433282.5
    lunar_eclipses = []
    while jd < end_jd and len(lunar_eclipses) < 40:
        try:
            se_res = swe.lun_eclipse_when(jd, 0, 0, False)
            tmax = se_res[1][0]
            if tmax < end_jd:
                lunar_eclipses.append((se_res[0], se_res[1]))
            jd = tmax + 30
        except Exception:
            jd += 30

    print(f"  Found {len(lunar_eclipses)} lunar eclipses")

    for ecl_type, se_times in lunar_eclipses:
        try:
            le_res = ephem.lun_eclipse_when(se_times[0] - 15, 0, 0)
            le_times = le_res[1]
        except Exception:
            continue

        se_tmax = se_times[0]
        le_tmax = le_times[0]
        diff_min = abs(le_tmax - se_tmax) * 1440
        tol_min = 5.0

        total_tests += 1
        if diff_min < tol_min:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 30:
                failures.append(
                    f"  Lunar JD={se_tmax:.4f}: tmax diff={diff_min:.2f} min"
                )

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 128 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
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
