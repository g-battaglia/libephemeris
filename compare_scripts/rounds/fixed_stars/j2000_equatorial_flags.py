#!/usr/bin/env python3
"""Round 161: Fixed stars with J2000+EQUATORIAL flags.

Compare fixed star positions with J2000, EQUATORIAL, and combined flag modes.
Tests the flag handling pipeline fixed in Round 95.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_J2000 = 32
FLG_EQUATORIAL = 2048
FLG_NONUT = 64
FLG_SIDEREAL = 65536

STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Altair",
    "Pollux",
    "Fomalhaut",
    "Deneb",
]

FLAG_COMBOS = [
    (FLG_SPEED, "default"),
    (FLG_SPEED | FLG_J2000, "J2000"),
    (FLG_SPEED | FLG_EQUATORIAL, "EQUATORIAL"),
    (FLG_SPEED | FLG_J2000 | FLG_EQUATORIAL, "J2000+EQUAT"),
    (FLG_SPEED | FLG_NONUT, "NONUT"),
    (FLG_SPEED | FLG_NONUT | FLG_EQUATORIAL, "NONUT+EQUAT"),
]

test_dates = []
for year in [1950, 1980, 2000, 2010, 2020, 2024, 2025, 2050]:
    jd = swe.julday(year, 1, 1, 12.0)
    test_dates.append((f"{year}", jd))

TOL_LON = 1.0
TOL_LAT = 1.0
TOL_SPD = 5.0  # Relaxed speed for known star speed diffs

passed = failed = errors = total = 0
failures = []

print(f"Round 161: Fixed Stars with J2000+EQUATORIAL Flags")
print(
    f"Testing {len(STARS)} stars x {len(FLAG_COMBOS)} flags x {len(test_dates)} dates"
)
print("=" * 90)

for label, jd in test_dates:
    for star_name in STARS:
        for flags, fname in FLAG_COMBOS:
            try:
                se_r = swe.fixstar2(star_name, jd, flags)
                se_d = se_r[0]

                le_r = ephem.fixstar2_ut(star_name, jd, flags)
                le_d = le_r[0]  # (pos_tuple, star_name, retflag)

                for i, (cn, mult, tol) in enumerate(
                    [
                        ("lon/RA", 3600, TOL_LON),
                        ("lat/Dec", 3600, TOL_LAT),
                        ("dist", 1, 100),
                        ("lon_spd", 3600, TOL_SPD),
                        ("lat_spd", 3600, TOL_SPD),
                        ("dist_spd", 1, 1e-3),
                    ]
                ):
                    total += 1
                    diff = abs(le_d[i] - se_d[i]) * mult
                    if i == 0 and diff > 180 * 3600:
                        diff = 360 * 3600 - diff
                    if diff <= tol:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f"  FAIL {label} {star_name} {fname} {cn}: SE={se_d[i]:.8f} LE={le_d[i]:.8f} diff={diff:.4f}"
                        )
                        if len(failures) <= 20:
                            print(failures[-1])
            except Exception as e:
                errors += 1
                if errors <= 3:
                    print(f"  ERROR {label} {star_name} {fname}: {e}")

print(f"\n{'=' * 90}")
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    cats = {}
    for f in failures:
        for _, fn in FLAG_COMBOS:
            if fn + ")" in f or fn + " " in f or f.endswith(fn):
                cats[fn] = cats.get(fn, 0) + 1
    print(f"\n{len(failures)} failures:")
    for c, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
else:
    print("\nAll tests passed!")
