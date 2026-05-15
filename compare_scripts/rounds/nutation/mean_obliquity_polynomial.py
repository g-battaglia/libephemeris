#!/usr/bin/env python3
"""Round 146: Mean obliquity polynomial precision at extreme dates.

Compare mean obliquity, true obliquity, and nutation components between
libephemeris and pyswisseph across a wide range of epochs, focusing on
extreme dates where polynomial divergence is most likely.

ECL_NUT (-1) returns: (true_obl, mean_obl, nut_lon, nut_obl, 0, 0)
"""

from __future__ import annotations
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

ECL_NUT = -1

# Test dates spanning from ~5000 BCE to ~3000 CE
# Focus on extreme dates + standard epochs
test_dates = []

# Extreme past (BCE dates)
for year in [-5000, -4000, -3000, -2500, -2000, -1500, -1000, -500]:
    jd = swe.julday(year, 1, 1, 12.0)
    test_dates.append((f"{year} BCE Jan 1", jd))

# Historical dates
for year in [0, 100, 500, 1000, 1200, 1400, 1500, 1600, 1700, 1800]:
    jd = swe.julday(year, 7, 1, 12.0)
    test_dates.append((f"{year} CE Jul 1", jd))

# Modern era (fine-grained)
for year in range(1900, 2030, 5):
    jd = swe.julday(year, 1, 1, 12.0)
    test_dates.append((f"{year} CE Jan 1", jd))

# J2000.0 exactly
test_dates.append(("J2000.0", 2451545.0))

# Future dates
for year in [2050, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]:
    jd = swe.julday(year, 1, 1, 12.0)
    test_dates.append((f"{year} CE Jan 1", jd))

# Extreme future
for year in [3500, 4000, 4500, 5000]:
    jd = swe.julday(year, 1, 1, 12.0)
    test_dates.append((f"{year} CE Jan 1", jd))

# Solstice/equinox dates for extra coverage
for year in [2000, 2024, 2025]:
    for month, day, label in [
        (3, 20, "VE"),
        (6, 21, "SS"),
        (9, 22, "AE"),
        (12, 21, "WS"),
    ]:
        jd = swe.julday(year, month, day, 12.0)
        test_dates.append((f"{year} {label}", jd))

# Tolerances in arcseconds
TOL_MEAN_OBL = 0.01  # 10 mas for mean obliquity
TOL_TRUE_OBL = 0.05  # 50 mas for true obliquity (includes nutation)
TOL_NUT_LON = 0.05  # 50 mas for nutation in longitude
TOL_NUT_OBL = 0.05  # 50 mas for nutation in obliquity
# Relax for extreme dates
TOL_EXTREME_MEAN = 1.0  # 1" for dates > 2000 years from J2000
TOL_EXTREME_NUT = 0.5  # 0.5" for nutation at extreme dates

passed = 0
failed = 0
errors = 0
total = 0
failures = []

print(f"Round 146: Mean Obliquity Polynomial Precision")
print(f"Testing {len(test_dates)} dates from ~5000 BCE to 5000 CE")
print("=" * 90)

for label, jd in test_dates:
    try:
        # pyswisseph: returns ((true_obl, mean_obl, nut_lon, nut_obl, 0, 0), retflag)
        se_result = swe.calc_ut(jd, ECL_NUT)
        se_data = se_result[0]
        se_true_obl = se_data[0]
        se_mean_obl = se_data[1]
        se_nut_lon = se_data[2]
        se_nut_obl = se_data[3]

        # libephemeris: same format
        le_result = ephem.calc_ut(jd, ECL_NUT, 0)
        le_data = le_result[0]
        le_true_obl = le_data[0]
        le_mean_obl = le_data[1]
        le_nut_lon = le_data[2]
        le_nut_obl = le_data[3]

        # Determine if extreme date (>2000 years from J2000)
        t_centuries = abs(jd - 2451545.0) / 36525.0
        is_extreme = t_centuries > 20.0  # >2000 years

        tol_mean = TOL_EXTREME_MEAN if is_extreme else TOL_MEAN_OBL
        tol_nut = TOL_EXTREME_NUT if is_extreme else TOL_NUT_LON
        tol_true = TOL_EXTREME_MEAN if is_extreme else TOL_TRUE_OBL

        # Compare mean obliquity (degrees -> arcseconds for diff)
        diff_mean = abs(le_mean_obl - se_mean_obl) * 3600.0
        diff_true = abs(le_true_obl - se_true_obl) * 3600.0
        diff_nut_lon = abs(le_nut_lon - se_nut_lon) * 3600.0
        diff_nut_obl = abs(le_nut_obl - se_nut_obl) * 3600.0

        tests = [
            ("mean_obl", diff_mean, tol_mean, se_mean_obl, le_mean_obl),
            ("true_obl", diff_true, tol_true, se_true_obl, le_true_obl),
            ("nut_lon", diff_nut_lon, tol_nut, se_nut_lon, le_nut_lon),
            ("nut_obl", diff_nut_obl, tol_nut, se_nut_obl, le_nut_obl),
        ]

        for name, diff, tol, se_val, le_val in tests:
            total += 1
            if diff <= tol:
                passed += 1
            else:
                failed += 1
                msg = f'  FAIL {label} {name}: SE={se_val:.10f} LE={le_val:.10f} diff={diff:.6f}" tol={tol}"'
                failures.append(msg)
                print(msg)

    except Exception as e:
        errors += 1
        total += 4
        print(f"  ERROR {label}: {e}")

print()
print("=" * 90)
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)

if failures:
    print(f"\nAll {len(failures)} failures:")
    for f in failures:
        print(f)
else:
    print("\nAll tests passed!")
