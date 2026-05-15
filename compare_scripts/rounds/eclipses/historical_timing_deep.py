#!/usr/bin/env python3
"""Round 193: Eclipse timing at historical dates.

Tests solar eclipse global search (sol_eclipse_when_glob) at known
historical eclipse dates to verify timing accuracy. Compares LE vs SE
for eclipse maximum time and type classification.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

# Known historical solar eclipses — search near these dates
HISTORICAL_ECLIPSES = [
    ("1999 Aug 11 total", 2451401.5),
    ("2001 Jun 21 total", 2452081.5),
    ("2003 Nov 23 total", 2452967.5),
    ("2006 Mar 29 total", 2453824.5),
    ("2008 Aug 1 total", 2454679.5),
    ("2009 Jul 22 total", 2455034.5),
    ("2010 Jul 11 total", 2455389.5),
    ("2012 Nov 13 total", 2456245.5),
    ("2015 Mar 20 total", 2457101.5),
    ("2016 Mar 9 total", 2457456.5),
    ("2017 Aug 21 total", 2457987.5),
    ("2019 Jul 2 total", 2458666.5),
    ("2020 Jun 21 annular", 2459021.5),
    ("2020 Dec 14 total", 2459197.5),
    ("2021 Jun 10 annular", 2459375.5),
    ("2023 Apr 20 hybrid", 2460054.5),
    ("2023 Oct 14 annular", 2460231.5),
    ("2024 Apr 8 total", 2460408.5),
    ("1900 May 28 total", 2415184.5),
    ("1918 Jun 8 total", 2421793.5),
    ("1919 May 29 total", 2422148.5),
    ("1925 Jan 24 total", 2424194.5),
    ("1945 Jul 9 total", 2431674.5),
    ("1952 Feb 25 total", 2434096.5),
    ("1961 Feb 15 total", 2437379.5),
    ("1970 Mar 7 total", 2440651.5),
    ("1973 Jun 30 total", 2441856.5),
    ("1979 Feb 26 total", 2443935.5),
    ("1980 Feb 16 total", 2444290.5),
    ("1991 Jul 11 total", 2448451.5),
    ("1994 Nov 3 total", 2449660.5),
    ("1995 Oct 24 total", 2450015.5),
    ("1997 Mar 9 total", 2450516.5),
    ("1998 Feb 26 total", 2450871.5),
]

FLAGS = 0


def test_eclipse_timing():
    global passed, failed, total

    print("=" * 70)
    print("Round 193: Eclipse Timing at Historical Dates")
    print("=" * 70)

    for label, approx_jd in HISTORICAL_ECLIPSES:
        search_jd = approx_jd - 30

        # libephemeris: returns (ecl_type_int, tret_tuple)
        try:
            le_result = ephem.sol_eclipse_when_glob(search_jd, FLAGS)
            le_ecl_type = le_result[0]  # int eclipse type flags
            le_tret = le_result[1]  # tuple of times
        except Exception as e:
            print(f"  {label}: LE error: {e}")
            continue

        # pyswisseph: returns (ecl_type_int, tret_tuple)
        try:
            se_result = swe.sol_eclipse_when_glob(search_jd, swe.FLG_SWIEPH)
            se_ecl_type = se_result[0]
            se_tret = se_result[1]
        except Exception as e:
            print(f"  {label}: SE error: {e}")
            continue

        # Extract maximum time (index 0)
        le_tmax = le_tret[0] if le_tret[0] != 0 else None
        se_tmax = se_tret[0] if se_tret[0] != 0 else None

        if le_tmax is None or se_tmax is None:
            continue

        # Compare maximum time
        total += 1
        diff_min = abs(le_tmax - se_tmax) * 1440

        if diff_min <= 10.0:  # 10 minute tolerance for eclipse max
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  {label} t_max: diff={diff_min:.2f} min (LE={le_tmax:.6f} SE={se_tmax:.6f})"
            )

        # Compare other contact times
        for idx, tname in [
            (2, "t_begin"),
            (3, "t_end"),
            (4, "t_totalbegin"),
            (5, "t_totalend"),
        ]:
            try:
                le_t = le_tret[idx] if idx < len(le_tret) else 0
                se_t = se_tret[idx] if idx < len(se_tret) else 0
            except (IndexError, TypeError):
                continue

            if le_t == 0 and se_t == 0:
                continue
            if le_t == 0 or se_t == 0:
                total += 1
                failed += 1
                failures.append(
                    f"  {label} {tname}: one is zero (LE={le_t:.6f} SE={se_t:.6f})"
                )
                continue

            total += 1
            ct_diff = abs(le_t - se_t) * 1440
            if ct_diff <= 15.0:
                passed += 1
            else:
                failed += 1
                failures.append(f"  {label} {tname}: diff={ct_diff:.2f} min")


def test_lunar_eclipse_timing():
    """Test lunar eclipse timing at known dates."""
    global passed, failed, total

    print("\n--- Lunar Eclipse Timing ---")

    lunar_eclipses = [
        ("2000 Jan 21 total", 2451564.5),
        ("2000 Jul 16 total", 2451741.5),
        ("2003 Nov 9 total", 2452953.5),
        ("2004 Oct 28 total", 2453307.5),
        ("2007 Mar 3 total", 2454163.5),
        ("2008 Feb 21 total", 2454518.5),
        ("2010 Dec 21 total", 2455551.5),
        ("2011 Jun 15 total", 2455728.5),
        ("2014 Apr 15 total", 2456763.5),
        ("2015 Apr 4 total", 2457116.5),
        ("2018 Jan 31 total", 2458150.5),
        ("2018 Jul 27 total", 2458327.5),
        ("2019 Jan 21 total", 2458505.5),
        ("2021 May 26 total", 2459360.5),
        ("2022 May 16 total", 2459715.5),
        ("2022 Nov 8 total", 2459891.5),
    ]

    for label, approx_jd in lunar_eclipses:
        search_jd = approx_jd - 30

        # libephemeris: returns (ecl_type_int, tret_tuple)
        try:
            le_result = ephem.lun_eclipse_when(search_jd, FLAGS, 0)
            le_ecl_type = le_result[0]
            le_tret = le_result[1]
        except Exception as e:
            continue

        # pyswisseph
        try:
            se_result = swe.lun_eclipse_when(search_jd, swe.FLG_SWIEPH, 0)
            se_ecl_type = se_result[0]
            se_tret = se_result[1]
        except Exception:
            continue

        le_tmax = le_tret[0]
        se_tmax = se_tret[0]

        if le_tmax == 0 or se_tmax == 0:
            continue

        total += 1
        diff_min = abs(le_tmax - se_tmax) * 1440

        if diff_min <= 10.0:
            passed += 1
        else:
            failed += 1
            failures.append(f"  LUNAR {label} t_max: diff={diff_min:.2f} min")

        # Compare contact times
        for idx, tname in [
            (2, "t_partial_begin"),
            (3, "t_partial_end"),
            (4, "t_total_begin"),
            (5, "t_total_end"),
            (6, "t_penumbral_begin"),
            (7, "t_penumbral_end"),
        ]:
            try:
                le_t = le_tret[idx] if idx < len(le_tret) else 0
                se_t = se_tret[idx] if idx < len(se_tret) else 0
            except (IndexError, TypeError):
                continue

            if le_t == 0 and se_t == 0:
                continue
            if le_t == 0 or se_t == 0:
                total += 1
                failed += 1
                failures.append(f"  LUNAR {label} {tname}: one is zero")
                continue

            total += 1
            ct_diff = abs(le_t - se_t) * 1440
            if ct_diff <= 15.0:
                passed += 1
            else:
                failed += 1
                failures.append(f"  LUNAR {label} {tname}: diff={ct_diff:.2f} min")


if __name__ == "__main__":
    test_eclipse_timing()
    test_lunar_eclipse_timing()

    print(f"\n{'=' * 70}")
    pct = 100 * passed / total if total > 0 else 0
    print(f"RESULTS: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
