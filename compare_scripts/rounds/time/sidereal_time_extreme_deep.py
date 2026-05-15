#!/usr/bin/env python3
"""Round 191: Sidereal time at extreme dates.

Tests sidereal time computation across a wide date range including
ancient and far-future dates, comparing against pyswisseph.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

# Test dates spanning a wide range
TEST_JDS = [
    # Ancient dates
    ("3000 BCE", 625673.5),
    ("2000 BCE", 990558.5),
    ("1000 BCE", 1355808.5),
    ("500 BCE", 1538395.5),
    ("0 CE (1 BCE)", 1721057.5),
    # Classical era
    ("500 CE", 1903682.5),
    ("1000 CE", 2086307.5),
    ("1500 CE", 2268932.5),
    # Modern era (fine sampling)
    ("1800 Jan 1", 2378497.0),
    ("1850 Jan 1", 2396758.5),
    ("1900 Jan 0.5", 2415020.0),
    ("1950 Jan 1", 2433282.5),
    ("1980 Jan 1", 2444239.5),
    ("1990 Jan 1", 2447892.5),
    ("2000 Jan 1.5 (J2000)", 2451545.0),
    ("2000 Jun 21", 2451716.5),
    ("2005 Jan 1", 2453371.5),
    ("2010 Jan 1", 2455197.5),
    ("2015 Jan 1", 2457023.5),
    ("2020 Jan 1", 2458849.5),
    ("2024 Mar 15", 2460384.5),
    # Future dates
    ("2050 Jan 1", 2469807.5),
    ("2100 Jan 1", 2488069.5),
    ("2200 Jan 1", 2524594.5),
    ("2500 Jan 1", 2634166.5),
    ("3000 Jan 1", 2816788.5),
]


def test_sidtime():
    global passed, failed, total

    print("=" * 70)
    print("Round 191: Sidereal Time at Extreme Dates")
    print("=" * 70)

    max_diff = 0.0
    max_diff_label = ""

    for label, jd in TEST_JDS:
        try:
            le_st = ephem.sidtime(jd)
            se_st = swe.sidtime(jd)
        except Exception as e:
            continue

        total += 1
        diff = abs(le_st - se_st)
        # Wrap around 24h
        if diff > 12:
            diff = 24 - diff
        diff_sec = diff * 3600  # hours to seconds

        if diff_sec > max_diff:
            max_diff = diff_sec
            max_diff_label = label

        # Tolerance: 0.01s for modern, 0.1s for ancient
        is_modern = jd > 2400000
        tol = 0.01 if is_modern else 0.1
        if diff_sec <= tol:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  {label} JD={jd}: LE={le_st:.10f}h SE={se_st:.10f}h diff={diff_sec:.4f}s"
            )

    print(f"\n  Max sidereal time diff: {max_diff:.6f}s at {max_diff_label}")

    # Also test sidtime0 (with obliquity and nutation)
    print("\n--- sidtime0 tests ---")
    for label, jd in TEST_JDS:
        if jd < 2300000:  # Skip very ancient for sidtime0
            continue
        try:
            # Get obliquity and nutation
            ecl_nut = ephem.calc_ut(jd, -1, 0)
            eps = ecl_nut[0][0]  # true obliquity
            nut_lon = ecl_nut[0][2]  # nutation in longitude

            le_st0 = ephem.sidtime0(jd, eps, nut_lon)
            se_st0 = swe.sidtime0(jd, eps, nut_lon)
        except Exception:
            continue

        total += 1
        diff = abs(le_st0 - se_st0)
        if diff > 12:
            diff = 24 - diff
        diff_sec = diff * 3600

        tol = 0.01
        if diff_sec <= tol:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  sidtime0 {label}: LE={le_st0:.10f}h SE={se_st0:.10f}h diff={diff_sec:.4f}s"
            )


def test_sidtime_at_greenwich_noon():
    """Test sidereal time at Greenwich noon for various dates."""
    global passed, failed, total

    print("\n--- Sidereal Time at Greenwich Noon ---")

    # Dense modern sampling
    for year_offset in range(0, 50):
        jd = 2451545.0 + year_offset * 365.25  # approximate yearly
        try:
            le_st = ephem.sidtime(jd)
            se_st = swe.sidtime(jd)
        except Exception:
            continue

        total += 1
        diff = abs(le_st - se_st)
        if diff > 12:
            diff = 24 - diff
        diff_sec = diff * 3600

        if diff_sec <= 0.005:
            passed += 1
        else:
            failed += 1
            failures.append(f"  noon JD={jd:.1f}: diff={diff_sec:.6f}s")


if __name__ == "__main__":
    test_sidtime()
    test_sidtime_at_greenwich_noon()

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
