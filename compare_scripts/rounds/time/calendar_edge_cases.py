#!/usr/bin/env python3
"""Round 199: Calendar edge cases (leap years, BCE dates, Julian/Gregorian boundary).

Tests julday and revjul at calendar edge cases including leap years,
century boundaries, BCE dates, and the Julian-Gregorian transition.
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

GREG_CAL = 1
JUL_CAL = 0


def test_julday():
    global passed, failed, total

    print("=" * 70)
    print("Round 199: Calendar Edge Cases")
    print("=" * 70)

    # Leap year dates
    print("\n--- Leap Year Dates ---")
    leap_dates = [
        (2000, 2, 29, 12.0, GREG_CAL, "2000 Feb 29 (leap)"),
        (2004, 2, 29, 12.0, GREG_CAL, "2004 Feb 29 (leap)"),
        (1900, 2, 28, 12.0, GREG_CAL, "1900 Feb 28 (century, no leap)"),
        (1600, 2, 29, 12.0, GREG_CAL, "1600 Feb 29 (400yr leap)"),
        (2024, 2, 29, 12.0, GREG_CAL, "2024 Feb 29 (leap)"),
        (2100, 2, 28, 12.0, GREG_CAL, "2100 Feb 28 (century, no leap)"),
        (2400, 2, 29, 12.0, GREG_CAL, "2400 Feb 29 (400yr leap)"),
    ]

    for y, m, d, h, cal, label in leap_dates:
        total += 1
        le_jd = ephem.julday(y, m, d, h, cal)
        se_jd = swe.julday(y, m, d, h, cal)
        diff = abs(le_jd - se_jd)
        if diff < 1e-10:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  julday {label}: LE={le_jd:.6f} SE={se_jd:.6f} diff={diff}"
            )

    # BCE dates
    print("\n--- BCE Dates ---")
    bce_dates = [
        (0, 1, 1, 12.0, JUL_CAL, "0 CE Jan 1 (=1 BCE) Julian"),
        (-1, 1, 1, 12.0, JUL_CAL, "-1 (=2 BCE) Jan 1 Julian"),
        (-4, 1, 1, 12.0, JUL_CAL, "-4 (=5 BCE) Jan 1 Julian"),
        (-100, 1, 1, 12.0, JUL_CAL, "-100 Jan 1 Julian"),
        (-500, 6, 15, 12.0, JUL_CAL, "-500 Jun 15 Julian"),
        (-1000, 1, 1, 12.0, JUL_CAL, "-1000 Jan 1 Julian"),
        (-3000, 1, 1, 12.0, JUL_CAL, "-3000 Jan 1 Julian"),
        (-4713, 1, 1, 12.0, JUL_CAL, "-4713 Jan 1 Julian (JD epoch)"),
        # Proleptic Gregorian
        (0, 1, 1, 12.0, GREG_CAL, "0 CE Jan 1 Gregorian"),
        (-1, 1, 1, 12.0, GREG_CAL, "-1 Jan 1 Gregorian"),
        (-100, 7, 4, 12.0, GREG_CAL, "-100 Jul 4 Gregorian"),
        (-1000, 3, 21, 12.0, GREG_CAL, "-1000 Mar 21 Gregorian"),
    ]

    for y, m, d, h, cal, label in bce_dates:
        total += 1
        le_jd = ephem.julday(y, m, d, h, cal)
        se_jd = swe.julday(y, m, d, h, cal)
        diff = abs(le_jd - se_jd)
        if diff < 1e-10:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  julday {label}: LE={le_jd:.6f} SE={se_jd:.6f} diff={diff}"
            )

    # Julian-Gregorian transition (Oct 1582)
    print("\n--- Julian-Gregorian Transition ---")
    transition_dates = [
        (1582, 10, 4, 12.0, JUL_CAL, "1582 Oct 4 Julian (last Julian day)"),
        (
            1582,
            10,
            15,
            12.0,
            GREG_CAL,
            "1582 Oct 15 Gregorian (first Gregorian day)",
        ),
        (1582, 10, 14, 12.0, JUL_CAL, "1582 Oct 14 Julian"),
        (1582, 10, 16, 12.0, GREG_CAL, "1582 Oct 16 Gregorian"),
        (1582, 1, 1, 12.0, JUL_CAL, "1582 Jan 1 Julian"),
        (1582, 12, 31, 12.0, GREG_CAL, "1582 Dec 31 Gregorian"),
    ]

    for y, m, d, h, cal, label in transition_dates:
        total += 1
        le_jd = ephem.julday(y, m, d, h, cal)
        se_jd = swe.julday(y, m, d, h, cal)
        diff = abs(le_jd - se_jd)
        if diff < 1e-10:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  julday {label}: LE={le_jd:.6f} SE={se_jd:.6f} diff={diff}"
            )


def test_revjul():
    global passed, failed, total

    print("\n--- revjul Roundtrip ---")

    # Test roundtrip: julday -> revjul -> julday
    test_jds = [
        0.0,
        0.5,
        1.0,
        1721057.5,  # 0 CE
        2299160.5,  # Just before Gregorian
        2299161.5,  # Just at Gregorian boundary
        2415020.0,  # 1900
        2451545.0,  # J2000
        2460310.5,  # 2024
        625673.5,  # 3000 BCE
        990558.5,  # 2000 BCE
    ]

    for jd in test_jds:
        for cal in [GREG_CAL, JUL_CAL]:
            cal_name = "Greg" if cal == GREG_CAL else "Jul"
            total += 1

            le_rev = ephem.revjul(jd, cal)
            se_rev = swe.revjul(jd, cal)

            # Compare year, month, day, hour
            le_y, le_m, le_d, le_h = le_rev
            se_y, se_m, se_d, se_h = se_rev

            if (
                le_y == se_y
                and le_m == se_m
                and le_d == se_d
                and abs(le_h - se_h) < 1e-10
            ):
                passed += 1
            else:
                failed += 1
                failures.append(
                    f"  revjul JD={jd} {cal_name}: LE=({le_y},{le_m},{le_d},{le_h:.6f}) SE=({se_y},{se_m},{se_d},{se_h:.6f})"
                )


def test_deltat_consistency():
    global passed, failed, total

    print("\n--- Delta-T Consistency ---")

    test_jds = [
        2451545.0,  # J2000
        2455197.5,  # 2010
        2458849.5,  # 2020
        2460310.5,  # 2024
        2415020.0,  # 1900
        2440000.0,  # 1968
    ]

    for jd in test_jds:
        total += 1
        le_dt = ephem.deltat(jd)
        se_dt = swe.deltat(jd)
        diff = abs(le_dt - se_dt) * 86400  # seconds

        if diff <= 1.0:  # 1 second tolerance
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  deltat JD={jd:.1f}: LE={le_dt:.10f} SE={se_dt:.10f} diff={diff:.4f}s"
            )


if __name__ == "__main__":
    test_julday()
    test_revjul()
    test_deltat_consistency()

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
