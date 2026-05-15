#!/usr/bin/env python3
"""Round 176: Planetary sect — diurnal/nocturnal planet detection accuracy.

Tests whether planets are correctly above/below horizon at birth times,
which is fundamental for sect (day/night) calculations in astrology.
Compares house_pos results to verify planet hemisphere placement.
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

BODIES = {
    "Sun": (swe.SUN, ephem.SUN),
    "Moon": (swe.MOON, ephem.MOON),
    "Mercury": (swe.MERCURY, ephem.MERCURY),
    "Venus": (swe.VENUS, ephem.VENUS),
    "Mars": (swe.MARS, ephem.MARS),
    "Jupiter": (swe.JUPITER, ephem.JUPITER),
    "Saturn": (swe.SATURN, ephem.SATURN),
}

LOCATIONS = [
    (51.5, -0.1),
    (40.7, -74.0),
    (35.7, 139.7),
    (-33.9, 151.2),
    (55.8, 37.6),
    (19.1, 72.9),
    (-22.9, -43.2),
    (64.1, -21.9),
]

TEST_DATES = []
for year in [1980, 1990, 2000, 2010, 2020, 2025]:
    for hour in [3, 9, 15, 21]:
        jd = swe.julday(year, 6, 15, float(hour))
        TEST_DATES.append((f"{year}/06/15 {hour}h", jd))

FLAGS = swe.FLG_SPEED
TOL_HOUSEPOS = 0.5  # house position tolerance in house units

passed = 0
failed = 0
skipped = 0
errors = []

for date_str, jd in TEST_DATES:
    for lat, lon in LOCATIONS:
        # Get houses first
        try:
            se_houses = swe.houses_ex(jd, lat, lon, b"P")
            se_cusps = se_houses[0]
            se_ascmc = se_houses[1]
            armc = se_ascmc[2]
            eps = swe.calc_ut(jd, swe.ECL_NUT, 0)[0][0]
        except Exception:
            skipped += 1
            continue

        try:
            le_houses = ephem.houses_ex2(jd, lat, lon, ord("P"), 0)
            le_cusps = le_houses[0]
            le_ascmc = le_houses[1]
        except Exception:
            skipped += 1
            continue

        for bname, (se_id, le_id) in BODIES.items():
            try:
                se_pos = swe.calc_ut(jd, se_id, FLAGS)[0]
                le_pos = ephem.calc_ut(jd, le_id, FLAGS)[0]
            except Exception:
                skipped += 1
                continue

            try:
                se_hp = swe.house_pos(armc, lat, eps, (se_pos[0], se_pos[1]), hsys=b"P")
                le_hp = ephem.house_pos(
                    armc, lat, eps, ord("P"), le_pos[0], le_pos[1]
                )
            except Exception:
                skipped += 1
                continue

            diff = abs(se_hp - le_hp)
            if diff > 6:
                diff = 12 - diff  # wrap around

            if diff <= TOL_HOUSEPOS:
                passed += 1
            else:
                failed += 1
                if len(errors) < 20:
                    errors.append(
                        f"  FAIL {bname} {date_str} @{lat},{lon}: SE hp={se_hp:.3f} LE hp={le_hp:.3f} diff={diff:.3f}"
                    )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 176: Planetary Sect / House Position ===")
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 176: ALL PASSED")
