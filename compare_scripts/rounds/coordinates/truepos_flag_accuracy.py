#!/usr/bin/env python3
"""Round 158: TRUEPOS flag accuracy deep.

Compare positions with FLG_TRUEPOS (geometric position without light-time
correction) for all planets across multiple epochs.
"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_TRUEPOS = 16
FLG_NOABERR = 1024

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

FLAG_COMBOS = [
    (FLG_SPEED | FLG_TRUEPOS, "TRUEPOS"),
    (FLG_SPEED | FLG_TRUEPOS | FLG_NOABERR, "TRUEPOS+NOABERR"),
]

test_dates = []
for year in range(1900, 2101, 10):
    jd = swe.julday(year, 6, 15, 12.0)
    test_dates.append((f"{year}", jd))
for year in [2024, 2025]:
    for month in range(1, 13):
        jd = swe.julday(year, month, 15, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

TOL_LON = 1.0
TOL_LAT = 1.0
TOL_DIST = 1e-6
TOL_SPD = 1.0

passed = failed = errors = total = 0
failures = []

print(f"Round 158: TRUEPOS Flag Accuracy Deep")
print(
    f"Testing {len(BODIES)} bodies x {len(FLAG_COMBOS)} flags x {len(test_dates)} dates"
)
print("=" * 90)

for label, jd in test_dates:
    for body, bname in BODIES.items():
        for flags, fname in FLAG_COMBOS:
            try:
                se_r = swe.calc_ut(jd, body, flags)[0]
                le_r = ephem.calc_ut(jd, body, flags)[0]
                for i, (cn, mult, tol) in enumerate(
                    [
                        ("lon", 3600, TOL_LON),
                        ("lat", 3600, TOL_LAT),
                        ("dist", 1, TOL_DIST),
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
                        failures.append(
                            f"  FAIL {label} {bname} {fname} {cn}: SE={se_r[i]:.8f} LE={le_r[i]:.8f} diff={diff:.6f}"
                        )
                        if len(failures) <= 10:
                            print(failures[-1])
            except Exception as e:
                errors += 1

print(f"\n{'=' * 90}")
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\n{len(failures)} failures")
    cats = {}
    for f in failures:
        for bn in BODIES.values():
            if bn in f:
                cats[bn] = cats.get(bn, 0) + 1
    for c, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
else:
    print("\nAll tests passed!")
