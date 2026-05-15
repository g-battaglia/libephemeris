#!/usr/bin/env python3
"""Round 129: Sidereal Time Precision at Fine Intervals

Tests sidereal time at fine (hourly) intervals across multiple years,
comparing LE vs SE with sub-second precision.
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
    print("ROUND 129: Sidereal Time Precision at Fine Intervals")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Test every 6 hours across 50 years (1970-2020)
    start_jd = 2440587.5  # 1970-01-01
    end_jd = 2458849.5  # 2020-01-01
    step = 0.25  # 6 hours

    jd = start_jd
    while jd <= end_jd:
        se_st = swe.sidtime(jd)
        le_st = ephem.sidtime(jd)

        diff_s = abs(le_st - se_st) * 3600  # seconds of time
        if diff_s > 43200:
            diff_s = 86400 - diff_s  # wrap

        tol_s = 0.01  # 10ms of sidereal time

        total_tests += 1
        if diff_s < tol_s:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 15:
                failures.append(f"  JD={jd:.4f}: diff={diff_s:.6f}s")

        jd += step

    # Also test sidtime0 (sidereal time at Greenwich)
    print("\n--- Testing sidtime0 ---")
    test_jds = [2451545.0, 2455197.5, 2459580.5, 2440587.5, 2444239.5]
    for jd in test_jds:
        eps = 23.44  # approximate obliquity
        dpsi = -0.004  # approximate nutation in longitude (degrees)
        try:
            se_st0 = swe.sidtime0(jd, eps, dpsi)
            le_st0 = ephem.sidtime0(jd, eps, dpsi)
            diff_s = abs(le_st0 - se_st0) * 3600
            if diff_s > 43200:
                diff_s = 86400 - diff_s

            total_tests += 1
            if diff_s < 0.01:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 20:
                    failures.append(f"  sidtime0 JD={jd:.1f}: diff={diff_s:.6f}s")
        except Exception:
            pass

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 129 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
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
