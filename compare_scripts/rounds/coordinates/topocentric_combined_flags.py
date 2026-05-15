#!/usr/bin/env python3
"""Round 160: Topocentric with NOABERR/NONUT combined flags.

Compare topocentric positions with various flag combinations to test
that parallax correction works correctly when nutation and/or aberration
are suppressed.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_TOPOCTR = 32768
FLG_NONUT = 64
FLG_NOABERR = 1024

BODIES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
}

FLAG_COMBOS = [
    (FLG_SPEED | FLG_TOPOCTR, "TOPO"),
    (FLG_SPEED | FLG_TOPOCTR | FLG_NONUT, "TOPO+NONUT"),
    (FLG_SPEED | FLG_TOPOCTR | FLG_NOABERR, "TOPO+NOABERR"),
    (FLG_SPEED | FLG_TOPOCTR | FLG_NONUT | FLG_NOABERR, "TOPO+NONUT+NOABERR"),
]

LOCATIONS = [
    (41.9, 12.5, 50, "Rome"),
    (51.5, -0.1, 11, "London"),
    (35.7, 139.7, 40, "Tokyo"),
    (-33.9, 151.2, 58, "Sydney"),
    (0.0, 0.0, 0, "Null_Island"),
]

test_dates = []
for year in [2000, 2010, 2020, 2024, 2025]:
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 1, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

TOL_LON = 1.0
TOL_LAT = 1.0
TOL_SPD = 5.0  # Topo Moon speed has known ~2-6"/day diff

passed = failed = errors = total = 0
failures = []

print(f"Round 160: Topocentric with NOABERR/NONUT Combined Flags")
print(
    f"Testing {len(BODIES)} bodies x {len(FLAG_COMBOS)} flags x {len(LOCATIONS)} locs x {len(test_dates)} dates"
)
print("=" * 90)

for lat, lon, alt, loc_name in LOCATIONS:
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)

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
                            ("dist_spd", 1, 1e-5),
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
                                f"  FAIL {loc_name} {label} {bname} {fname} {cn}: diff={diff:.4f}"
                            )
                            if len(failures) <= 15:
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
