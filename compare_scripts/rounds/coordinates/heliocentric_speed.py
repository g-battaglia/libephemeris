#!/usr/bin/env python3
"""Round 159: Heliocentric speed accuracy.

Compare heliocentric speeds for all planets. Heliocentric mode removes
geocentric parallax effects, so speed differences reveal pure ephemeris
model differences.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_HELCTR = 8
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
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
}

test_dates = []
for year in range(1950, 2051, 5):
    jd = swe.julday(year, 6, 15, 12.0)
    test_dates.append((f"{year}", jd))

TOL_LON = 1.0
TOL_LAT = 1.0
TOL_SPD = 1.0

passed = failed = errors = total = 0
failures = []

print(f"Round 159: Heliocentric Speed Accuracy")
print(f"Testing {len(BODIES)} bodies x {len(test_dates)} dates")
print("=" * 90)

for label, jd in test_dates:
    for body, bname in BODIES.items():
        if body in (0, 1):
            continue  # Sun/Moon not meaningful in helio
        try:
            se_r = swe.calc_ut(jd, body, FLG_SPEED | FLG_HELCTR)[0]
            le_r = ephem.calc_ut(jd, body, FLG_SPEED | FLG_HELCTR)[0]
            for i, (cn, mult, tol) in enumerate(
                [
                    ("lon", 3600, TOL_LON),
                    ("lat", 3600, TOL_LAT),
                    ("dist", 1, 1e-6),
                    ("lon_spd", 3600, TOL_SPD),
                    ("lat_spd", 3600, TOL_SPD),
                    ("dist_spd", 1, 1e-6),
                ]
            ):
                total += 1
                diff = abs(le_r[i] - se_r[i]) * mult
                if i == 0 and diff > 180 * 3600:
                    diff = 360 * 3600 - diff
                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    failures.append(f"  FAIL {label} {bname} {cn}: diff={diff:.6f}")
                    if len(failures) <= 10:
                        print(failures[-1])
        except Exception as e:
            errors += 1

print(f"\n{'=' * 90}")
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    cats = {}
    for f in failures:
        for bn in BODIES.values():
            if bn in f:
                cats[bn] = cats.get(bn, 0) + 1
    print(f"\n{len(failures)} failures:")
    for c, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
else:
    print("\nAll tests passed!")
