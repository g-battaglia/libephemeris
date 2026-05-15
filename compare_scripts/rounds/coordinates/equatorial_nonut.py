#!/usr/bin/env python3
"""Round 153: Equatorial RA/Dec with NONUT flag.

Compare equatorial coordinates (FLG_EQUATORIAL) with and without NONUT flag
for all planets across multiple epochs. NONUT removes nutation to give
mean equatorial coordinates instead of true equatorial.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_EQUATORIAL = 2048
FLG_NONUT = 64

BODIES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

# Flag combinations
FLAG_COMBOS = [
    (FLG_SPEED | FLG_EQUATORIAL, "EQUATORIAL"),
    (FLG_SPEED | FLG_EQUATORIAL | FLG_NONUT, "EQUATORIAL+NONUT"),
    (FLG_SPEED | FLG_NONUT, "NONUT_ecliptic"),
]

test_dates = []
for year in [1900, 1950, 1980, 2000, 2010, 2020, 2024, 2025, 2050, 2100]:
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 1, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

TOL_LON = 1.0  # arcsec for RA/lon
TOL_LAT = 1.0  # arcsec for Dec/lat
TOL_SPD = 1.0  # "/day

passed = failed = errors = total = 0
failures = []

print(f"Round 153: Equatorial RA/Dec with NONUT Flag")
print(
    f"Testing {len(BODIES)} bodies x {len(FLAG_COMBOS)} flags x {len(test_dates)} dates"
)
print("=" * 90)

for label, jd in test_dates:
    for body, bname in BODIES.items():
        for flags, fname in FLAG_COMBOS:
            try:
                se_r = swe.calc_ut(jd, body, flags)
                le_r = ephem.calc_ut(jd, body, flags)
                se_d = se_r[0]
                le_d = le_r[0]

                comp_names = [
                    "RA/lon",
                    "Dec/lat",
                    "dist",
                    "RA_spd",
                    "Dec_spd",
                    "dist_spd",
                ]
                for i in range(6):
                    total += 1
                    se_val = se_d[i]
                    le_val = le_d[i]

                    if i <= 1 or i == 3 or i == 4:
                        diff = abs(le_val - se_val) * 3600.0
                        # Handle RA wrap
                        if i == 0 and diff > 180 * 3600:
                            diff = 360 * 3600 - diff
                        tol = TOL_LON if i <= 1 else TOL_SPD
                    elif i == 2:
                        diff = abs(le_val - se_val)
                        tol = 1e-6
                    else:
                        diff = abs(le_val - se_val)
                        tol = 1e-6

                    if diff <= tol:
                        passed += 1
                    else:
                        failed += 1
                        msg = f"  FAIL {label} {bname} {fname} {comp_names[i]}: SE={se_val:.8f} LE={le_val:.8f} diff={diff:.6f}"
                        failures.append(msg)
                        if len(failures) <= 20:
                            print(msg)
            except Exception as e:
                errors += 1

print()
print("=" * 90)
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\n{len(failures)} failures total")
    cats = {}
    for f in failures:
        for _, fname in FLAG_COMBOS:
            if fname in f:
                cats[fname] = cats.get(fname, 0) + 1
    for c, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
else:
    print("\nAll tests passed!")
