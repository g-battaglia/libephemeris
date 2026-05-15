"""Round 138: Julian Day Conversion Round-Trip Integrity.

Tests that julday -> revjul -> julday round-trips perfectly for all
calendar types, including edge cases (BCE dates, Gregorian transition,
leap years, midnight vs noon, fractional days).
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

GREG_CAL = 1
JUL_CAL = 0

passed = 0
failed = 0
total = 0

# Test 1: julday -> revjul -> julday round-trip
print("=== Test 1: julday -> revjul -> julday round-trip ===")

test_dates = [
    # (year, month, day, hour, gregflag, label)
    (-4712, 1, 1, 12.0, JUL_CAL, "JD epoch"),
    (-3000, 6, 15, 6.5, JUL_CAL, "3001 BCE"),
    (-1000, 3, 21, 12.0, JUL_CAL, "1001 BCE vernal equinox"),
    (-500, 12, 31, 23.99, JUL_CAL, "501 BCE year end"),
    (-1, 1, 1, 0.0, JUL_CAL, "2 BCE"),
    (0, 1, 1, 0.0, JUL_CAL, "1 BCE astronomical"),
    (1, 1, 1, 0.0, JUL_CAL, "1 CE start"),
    (100, 7, 4, 18.0, JUL_CAL, "100 CE"),
    (325, 6, 20, 12.0, JUL_CAL, "Council of Nicaea"),
    (1000, 1, 1, 0.0, JUL_CAL, "1000 CE Julian"),
    (1582, 10, 4, 12.0, JUL_CAL, "Last Julian day"),
    (1582, 10, 15, 12.0, GREG_CAL, "First Gregorian day"),
    (1582, 10, 15, 0.0, GREG_CAL, "Gregorian midnight"),
    (1600, 2, 29, 12.0, GREG_CAL, "Gregorian leap day"),
    (1700, 2, 28, 12.0, GREG_CAL, "Not leap in Gregorian"),
    (1900, 1, 1, 0.0, GREG_CAL, "1900 start"),
    (1900, 2, 28, 12.0, GREG_CAL, "1900 not leap"),
    (1969, 12, 31, 23.999, GREG_CAL, "Pre-Unix epoch"),
    (1970, 1, 1, 0.0, GREG_CAL, "Unix epoch"),
    (2000, 1, 1, 12.0, GREG_CAL, "J2000"),
    (2000, 2, 29, 0.0, GREG_CAL, "Y2K leap day"),
    (2024, 3, 20, 3.06, GREG_CAL, "Vernal equinox 2024"),
    (2024, 6, 15, 12.0, GREG_CAL, "Current epoch"),
    (2100, 2, 28, 12.0, GREG_CAL, "2100 not leap"),
    (3000, 12, 31, 23.999, GREG_CAL, "Year 3000"),
    # Proleptic Gregorian
    (100, 1, 1, 12.0, GREG_CAL, "100 CE proleptic Greg"),
    (-500, 6, 15, 12.0, GREG_CAL, "501 BCE proleptic Greg"),
    (-4000, 1, 1, 0.0, GREG_CAL, "4001 BCE proleptic Greg"),
]

for year, month, day, hour, gregflag, label in test_dates:
    total += 1
    try:
        # SE round-trip
        se_jd = swe.julday(year, month, day, hour, gregflag)
        se_rev = swe.revjul(se_jd, gregflag)
        se_jd2 = swe.julday(se_rev[0], se_rev[1], se_rev[2], se_rev[3], gregflag)

        # LE round-trip
        le_jd = ephem.julday(year, month, day, hour, gregflag)
        le_rev = ephem.revjul(le_jd, gregflag)
        le_jd2 = ephem.julday(le_rev[0], le_rev[1], le_rev[2], le_rev[3], gregflag)

        # Check JD match
        jd_diff = abs(se_jd - le_jd)
        rt_diff_se = abs(se_jd - se_jd2)
        rt_diff_le = abs(le_jd - le_jd2)

        # JD should match exactly
        if jd_diff > 1e-10:
            failed += 1
            print(
                f"FAIL JD mismatch {label}: SE={se_jd:.10f} LE={le_jd:.10f} diff={jd_diff:.2e}"
            )
        elif rt_diff_le > 1e-8:
            failed += 1
            print(
                f"FAIL RT mismatch {label}: LE_JD={le_jd:.10f} LE_RT={le_jd2:.10f} diff={rt_diff_le:.2e}"
            )
        elif (
            abs(se_rev[0] - le_rev[0]) > 0
            or abs(se_rev[1] - le_rev[1]) > 0
            or abs(se_rev[2] - le_rev[2]) > 0
        ):
            failed += 1
            print(f"FAIL revjul date mismatch {label}: SE={se_rev[:3]} LE={le_rev[:3]}")
        elif abs(se_rev[3] - le_rev[3]) > 1e-6:
            failed += 1
            print(
                f"FAIL revjul hour mismatch {label}: SE_h={se_rev[3]:.8f} LE_h={le_rev[3]:.8f}"
            )
        else:
            passed += 1
    except Exception as e:
        failed += 1
        print(f"ERR  {label}: {e}")

# Test 2: revjul -> julday -> revjul round-trip from JD values
print("\n=== Test 2: JD -> revjul -> julday -> revjul consistency ===")

import random

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

random.seed(42)

# Generate random JDs spanning full range
test_jds = [
    0.0,  # 4713 BCE
    100000.0,  # ~4440 BCE
    500000.0,  # ~3344 BCE
    1000000.0,  # ~1974 BCE
    1500000.0,  # ~604 BCE
    1721057.5,  # 1 BCE Jan 1 Julian
    1721425.5,  # 1 CE Jan 1
    2000000.0,  # ~763 CE
    2299160.5,  # Oct 4, 1582 (last Julian)
    2299161.5,  # Oct 15, 1582 (first Gregorian)
    2415020.5,  # 1900 Jan 1
    2440587.5,  # 1970 Jan 1
    2451545.0,  # 2000 Jan 1.5 (J2000)
    2460676.5,  # 2025 Jan 1
    2488069.5,  # 2100 Jan 1
    2816787.5,  # 3000 Jan 1
]

# Add random fractional JDs
for _ in range(50):
    test_jds.append(random.uniform(0.0, 3000000.0))

for jd in test_jds:
    for gregflag in [GREG_CAL, JUL_CAL]:
        total += 1
        try:
            se_rev = swe.revjul(jd, gregflag)
            se_jd2 = swe.julday(se_rev[0], se_rev[1], se_rev[2], se_rev[3], gregflag)

            le_rev = ephem.revjul(jd, gregflag)
            le_jd2 = ephem.julday(
                le_rev[0], le_rev[1], le_rev[2], le_rev[3], gregflag
            )

            # Revjul should give same date components
            date_match = (
                se_rev[0] == le_rev[0]
                and se_rev[1] == le_rev[1]
                and se_rev[2] == le_rev[2]
            )
            hour_match = abs(se_rev[3] - le_rev[3]) < 1e-6

            # Round-trip JD should be close to original
            se_rt_err = abs(jd - se_jd2)
            le_rt_err = abs(jd - le_jd2)

            if not date_match:
                failed += 1
                cal = "Greg" if gregflag else "Jul"
                print(f"FAIL date {cal} JD={jd:.1f}: SE={se_rev[:3]} LE={le_rev[:3]}")
            elif not hour_match:
                failed += 1
                cal = "Greg" if gregflag else "Jul"
                print(
                    f"FAIL hour {cal} JD={jd:.1f}: SE_h={se_rev[3]:.8f} LE_h={le_rev[3]:.8f}"
                )
            elif le_rt_err > 1e-8:
                failed += 1
                cal = "Greg" if gregflag else "Jul"
                print(
                    f"FAIL RT {cal} JD={jd:.1f}: orig={jd:.10f} rt={le_jd2:.10f} err={le_rt_err:.2e}"
                )
            else:
                passed += 1
        except Exception as e:
            failed += 1
            cal = "Greg" if gregflag else "Jul"
            print(f"ERR  {cal} JD={jd:.1f}: {e}")

# Test 3: Delta-T consistency
print("\n=== Test 3: Delta-T consistency ===")

deltat_jds = [
    2415020.5,  # 1900
    2433282.5,  # 1950
    2440587.5,  # 1970
    2444239.5,  # 1980
    2447892.5,  # 1990
    2451545.0,  # 2000
    2455197.5,  # 2010
    2458849.5,  # 2020
    2460676.5,  # 2025
]

for jd in deltat_jds:
    total += 1
    try:
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff_sec = abs(se_dt - le_dt) * 86400.0  # convert days to seconds

        if diff_sec < 1.0:  # 1 second tolerance
            passed += 1
        else:
            failed += 1
            print(
                f"FAIL DeltaT JD={jd}: SE={se_dt * 86400:.3f}s LE={le_dt * 86400:.3f}s diff={diff_sec:.3f}s"
            )
    except Exception as e:
        failed += 1
        print(f"ERR  DeltaT JD={jd}: {e}")

print(f"\n{'=' * 60}")
print(f"Round 138: JD Conversion Round-Trip Integrity")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
