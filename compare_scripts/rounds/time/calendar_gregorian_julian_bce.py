#!/usr/bin/env python3
"""Round 104: Calendar Conversion Edge Cases

Tests julday/revjul at calendar boundaries and edge cases:
1. Gregorian/Julian calendar transition (Oct 15, 1582)
2. BCE dates (negative years)
3. Leap year boundaries
4. Proleptic Gregorian calendar
5. Round-trip consistency (date -> JD -> date)
6. Extreme dates (ancient and far future)
"""

from __future__ import annotations
import sys, os, time
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

GREG_CAL = 1
JUL_CAL = 0


def run_tests():
    passed = 0
    failed = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 104: Calendar Conversion Edge Cases")
    print("=" * 80)

    # =========================================================================
    # PART 1: julday agreement — common dates
    # =========================================================================
    print("\n--- PART 1: julday Common Dates ---")
    p1_pass = 0
    p1_fail = 0

    common_dates = [
        (2000, 1, 1, 12.0, GREG_CAL, "J2000"),
        (2000, 1, 1, 0.0, GREG_CAL, "2000-01-01 0h"),
        (1900, 1, 1, 12.0, GREG_CAL, "1900-01-01"),
        (1970, 1, 1, 0.0, GREG_CAL, "Unix epoch"),
        (2024, 3, 20, 12.0, GREG_CAL, "2024 equinox"),
        (1582, 10, 15, 0.0, GREG_CAL, "Greg start"),
        (1582, 10, 4, 0.0, JUL_CAL, "Julian end"),
        (100, 1, 1, 0.0, JUL_CAL, "100 AD Julian"),
        (-4712, 1, 1, 12.0, JUL_CAL, "JD epoch"),
        (2100, 12, 31, 23.999, GREG_CAL, "Far future"),
    ]

    for year, month, day, hour, cal, label in common_dates:
        total += 1
        se_jd = swe.julday(year, month, day, hour, cal)
        le_jd = ephem.julday(year, month, day, hour, cal)
        diff = abs(se_jd - le_jd)

        if diff < 1e-10:
            p1_pass += 1
            passed += 1
        else:
            p1_fail += 1
            failed += 1
            if len(fail_details) < 10:
                fail_details.append(
                    f"  FAIL [julday] {label}: SE={se_jd:.10f} LE={le_jd:.10f} diff={diff:.2e}"
                )

    print(f"  Part 1: {p1_pass}/{p1_pass + p1_fail} passed")

    # =========================================================================
    # PART 2: revjul agreement — JD -> date
    # =========================================================================
    print("\n--- PART 2: revjul JD -> Date ---")
    p2_pass = 0
    p2_fail = 0

    test_jds = [
        (2451545.0, GREG_CAL, "J2000"),
        (2440000.0, GREG_CAL, "1968"),
        (2415020.0, GREG_CAL, "1900"),
        (2299161.0, GREG_CAL, "Greg start JD"),
        (2299160.0, JUL_CAL, "Julian end JD"),
        (0.0, JUL_CAL, "JD 0"),
        (1000000.0, JUL_CAL, "JD 1M"),
        (2500000.0, GREG_CAL, "JD 2.5M"),
        (1721425.5, JUL_CAL, "1 AD Jan 1"),
        (2460000.0, GREG_CAL, "2023"),
    ]

    for jd, cal, label in test_jds:
        total += 1
        se_result = swe.revjul(jd, cal)
        le_result = ephem.revjul(jd, cal)

        # SE returns (year, month, day, hour)
        # LE returns (year, month, day, hour)
        se_y, se_m, se_d, se_h = se_result
        le_y, le_m, le_d, le_h = le_result

        if se_y == le_y and se_m == le_m and se_d == le_d and abs(se_h - le_h) < 1e-8:
            p2_pass += 1
            passed += 1
        else:
            p2_fail += 1
            failed += 1
            if len(fail_details) < 15:
                fail_details.append(
                    f"  FAIL [revjul] {label} JD={jd}: "
                    f"SE=({se_y},{se_m},{se_d},{se_h:.6f}) "
                    f"LE=({le_y},{le_m},{le_d},{le_h:.6f})"
                )

    print(f"  Part 2: {p2_pass}/{p2_pass + p2_fail} passed")

    # =========================================================================
    # PART 3: Round-trip (date -> JD -> date)
    # =========================================================================
    print("\n--- PART 3: Round-Trip (date -> JD -> date) ---")
    p3_pass = 0
    p3_fail = 0

    roundtrip_dates = [
        (2024, 2, 29, 12.0, GREG_CAL, "Leap day 2024"),
        (1900, 2, 28, 12.0, GREG_CAL, "Non-leap 1900"),
        (2000, 2, 29, 12.0, GREG_CAL, "Leap day 2000"),
        (1, 1, 1, 0.0, JUL_CAL, "1 AD"),
        (-1, 1, 1, 0.0, JUL_CAL, "2 BC"),
        (-4712, 1, 1, 12.0, JUL_CAL, "JD epoch"),
        (3000, 6, 15, 6.5, GREG_CAL, "Far future"),
        (-3000, 7, 1, 12.0, JUL_CAL, "Ancient"),
    ]

    for year, month, day, hour, cal, label in roundtrip_dates:
        total += 1
        # LE round-trip
        le_jd = ephem.julday(year, month, day, hour, cal)
        le_y, le_m, le_d, le_h = ephem.revjul(le_jd, cal)

        if le_y == year and le_m == month and le_d == day and abs(le_h - hour) < 1e-8:
            p3_pass += 1
            passed += 1
        else:
            p3_fail += 1
            failed += 1
            if len(fail_details) < 20:
                fail_details.append(
                    f"  FAIL [RT] {label}: in=({year},{month},{day},{hour}) "
                    f"out=({le_y},{le_m},{le_d},{le_h:.6f})"
                )

    print(f"  Part 3: {p3_pass}/{p3_pass + p3_fail} passed")

    # =========================================================================
    # PART 4: Proleptic Gregorian calendar dates
    # =========================================================================
    print("\n--- PART 4: Proleptic Gregorian (pre-1582) ---")
    p4_pass = 0
    p4_fail = 0

    proleptic_dates = [
        (1000, 6, 15, 12.0, GREG_CAL, "1000 AD Greg"),
        (500, 1, 1, 12.0, GREG_CAL, "500 AD Greg"),
        (1, 1, 1, 12.0, GREG_CAL, "1 AD Greg"),
        (-500, 1, 1, 12.0, GREG_CAL, "-500 Greg"),
        (-1000, 6, 15, 12.0, GREG_CAL, "-1000 Greg"),
    ]

    for year, month, day, hour, cal, label in proleptic_dates:
        total += 1
        se_jd = swe.julday(year, month, day, hour, cal)
        le_jd = ephem.julday(year, month, day, hour, cal)
        diff = abs(se_jd - le_jd)

        if diff < 1e-10:
            p4_pass += 1
            passed += 1
        else:
            p4_fail += 1
            failed += 1
            if len(fail_details) < 25:
                fail_details.append(
                    f"  FAIL [proleptic] {label}: SE={se_jd:.6f} LE={le_jd:.6f} diff={diff:.2e}"
                )

    print(f"  Part 4: {p4_pass}/{p4_pass + p4_fail} passed")

    # =========================================================================
    # PART 5: BCE dates with negative years
    # =========================================================================
    print("\n--- PART 5: BCE Dates (negative years) ---")
    p5_pass = 0
    p5_fail = 0

    bce_dates = [
        (-4713, 11, 24, 12.0, JUL_CAL, "JD 0 date"),
        (-30, 6, 15, 12.0, JUL_CAL, "31 BC"),
        (-100, 3, 1, 0.0, JUL_CAL, "101 BC"),
        (-1000, 12, 31, 23.5, JUL_CAL, "1001 BC"),
        (-5000, 1, 1, 12.0, JUL_CAL, "5001 BC"),
    ]

    for year, month, day, hour, cal, label in bce_dates:
        total += 1
        se_jd = swe.julday(year, month, day, hour, cal)
        le_jd = ephem.julday(year, month, day, hour, cal)
        diff = abs(se_jd - le_jd)

        if diff < 1e-10:
            p5_pass += 1
            passed += 1
        else:
            p5_fail += 1
            failed += 1
            if len(fail_details) < 30:
                fail_details.append(
                    f"  FAIL [BCE] {label}: SE={se_jd:.6f} LE={le_jd:.6f} diff={diff:.2e}"
                )

        # Also test revjul round-trip
        total += 1
        le_y, le_m, le_d, le_h = ephem.revjul(le_jd, cal)
        se_y, se_m, se_d, se_h = swe.revjul(se_jd, cal)

        if se_y == le_y and se_m == le_m and se_d == le_d and abs(se_h - le_h) < 1e-6:
            p5_pass += 1
            passed += 1
        else:
            p5_fail += 1
            failed += 1
            if len(fail_details) < 35:
                fail_details.append(
                    f"  FAIL [BCE_RT] {label}: "
                    f"SE=({se_y},{se_m},{se_d},{se_h:.4f}) "
                    f"LE=({le_y},{le_m},{le_d},{le_h:.4f})"
                )

    print(f"  Part 5: {p5_pass}/{p5_pass + p5_fail} passed")

    # =========================================================================
    # PART 6: Delta-T at calendar boundaries
    # =========================================================================
    print("\n--- PART 6: Delta-T at Calendar Boundaries ---")
    p6_pass = 0
    p6_fail = 0

    deltat_jds = [
        2451545.0,  # J2000
        2440000.0,  # 1968
        2415020.0,  # 1900
        2460000.0,  # 2023
        2299161.0,  # Greg start
    ]

    for jd in deltat_jds:
        total += 1
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff = abs(se_dt - le_dt) * 86400  # seconds

        if diff < 1.0:  # 1 second tolerance
            p6_pass += 1
            passed += 1
        else:
            p6_fail += 1
            failed += 1
            if len(fail_details) < 40:
                fail_details.append(
                    f"  FAIL [deltat] JD={jd:.1f}: SE={se_dt:.10f} LE={le_dt:.10f} diff={diff:.3f}s"
                )

    print(f"  Part 6: {p6_pass}/{p6_pass + p6_fail} passed")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    pct = 100 * passed / max(1, total)
    print()
    print("=" * 80)
    print(f"ROUND 104 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
