#!/usr/bin/env python3
"""Round 157: Sidereal mode with NONUT flag.

Compare positions with FLG_SIDEREAL + FLG_NONUT combination for all planets.
This tests that sidereal ayanamsha subtraction works correctly when nutation
is suppressed (mean ecliptic of date).
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
FLG_NONUT = 64
FLG_SIDEREAL = 65536

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

# Test several sidereal modes
SIDM_MODES = [0, 1, 3, 27]  # Fagan/Bradley, Lahiri, Raman, True Citra
SIDM_NAMES = {0: "Fagan", 1: "Lahiri", 3: "Raman", 27: "TrueCitra"}

FLAG_COMBOS = [
    (FLG_SPEED | FLG_SIDEREAL, "SIDEREAL"),
    (FLG_SPEED | FLG_SIDEREAL | FLG_NONUT, "SIDEREAL+NONUT"),
]

test_dates = []
for year in [1950, 1980, 2000, 2010, 2020, 2024, 2025, 2050]:
    for month in [1, 7]:
        jd = swe.julday(year, month, 1, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

TOL_LON = 15.0  # 15" (includes ~14" sidereal model difference)
TOL_LAT = 1.0
TOL_SPD = 1.0

passed = failed = errors = total = 0
failures = []

print(f"Round 157: Sidereal Mode with NONUT Flag")
print(
    f"Testing {len(BODIES)} bodies x {len(SIDM_MODES)} modes x {len(FLAG_COMBOS)} flags x {len(test_dates)} dates"
)
print("=" * 90)

for sidm in SIDM_MODES:
    swe.set_sid_mode(sidm)
    ephem.set_sid_mode(sidm, 0, 0)

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
                            failures.append(
                                f"  FAIL {SIDM_NAMES[sidm]} {label} {bname} {fname} {cn}: diff={diff:.4f}"
                            )
                            if len(failures) <= 15:
                                print(failures[-1])
                except Exception as e:
                    errors += 1

# Reset sidereal mode
swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0, 0)

print(f"\n{'=' * 90}")
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    cats = {}
    for f in failures:
        for _, fn in FLAG_COMBOS:
            if fn in f:
                cats[fn] = cats.get(fn, 0) + 1
    print(f"\n{len(failures)} failures:")
    for c, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
else:
    print("\nAll tests passed!")
