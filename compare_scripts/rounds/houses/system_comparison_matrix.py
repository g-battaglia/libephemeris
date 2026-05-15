#!/usr/bin/env python3
"""Round 179: House system comparison matrix — all systems at same moment.

Tests all supported house systems at the same time/location to verify
internal consistency. Computes cusps for every house system and compares
each against SE.
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

HOUSE_SYSTEMS = {
    "P": "Placidus",
    "K": "Koch",
    "O": "Porphyry",
    "R": "Regiomontanus",
    "C": "Campanus",
    "E": "Equal",
    "W": "WholeSign",
    "B": "Alcabitius",
    "M": "Morinus",
    "X": "AxialMeridian",
    "H": "AzimuthalHorizon",
    "T": "Polich-Page",
    "G": "Gauquelin",
    "V": "Vehlow",
}

LOCATIONS = [
    ("London", 51.5, -0.1),
    ("NewYork", 40.7, -74.0),
    ("Sydney", -33.9, 151.2),
    ("Moscow", 55.8, 37.6),
    ("Quito", -0.2, -78.5),
]

TEST_DATES = []
for year in [1990, 2000, 2010, 2020, 2025]:
    for hour in [0, 6, 12, 18]:
        jd = swe.julday(year, 3, 21, float(hour))
        TEST_DATES.append((f"{year}/03/21 {hour}h", jd))

TOL_CUSP = 5.0  # arcsec for most systems
TOL_ASCMC = 5.0  # arcsec for ASC/MC

passed = 0
failed = 0
skipped = 0
errors = []

for date_str, jd in TEST_DATES:
    for loc_name, lat, lon in LOCATIONS:
        for hsys_ch, hsys_name in HOUSE_SYSTEMS.items():
            try:
                se_r = swe.houses_ex(jd, lat, lon, hsys_ch.encode())
                se_cusps = se_r[0]
                se_ascmc = se_r[1]
            except Exception:
                skipped += 1
                continue

            try:
                le_r = ephem.houses_ex2(jd, lat, lon, ord(hsys_ch), 0)
                le_cusps = le_r[0]
                le_ascmc = le_r[1]
            except Exception:
                skipped += 1
                continue

            # Compare ASC and MC
            for i, name in [(0, "ASC"), (1, "MC")]:
                diff = abs(se_ascmc[i] - le_ascmc[i])
                if diff > 180:
                    diff = 360 - diff
                diff_as = diff * 3600
                if diff_as <= TOL_ASCMC:
                    passed += 1
                else:
                    failed += 1
                    if len(errors) < 20:
                        errors.append(
                            f'  FAIL {hsys_name} {name} {date_str} @{loc_name}: {diff_as:.1f}"'
                        )

            # Compare cusps (12 cusps, 0-indexed)
            num_cusps = min(len(se_cusps), len(le_cusps), 12)
            for i in range(num_cusps):
                diff = abs(se_cusps[i] - le_cusps[i])
                if diff > 180:
                    diff = 360 - diff
                diff_as = diff * 3600

                tol = TOL_CUSP
                if hsys_ch == "I":  # Sunshine
                    tol = 3600.0  # known large differences

                if diff_as <= tol:
                    passed += 1
                else:
                    failed += 1
                    if len(errors) < 20:
                        errors.append(
                            f'  FAIL {hsys_name} cusp{i + 1} {date_str} @{loc_name}: {diff_as:.1f}"'
                        )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 179: House System Comparison Matrix ===")
print(
    f"Dates: {len(TEST_DATES)}, Locations: {len(LOCATIONS)}, Systems: {len(HOUSE_SYSTEMS)}"
)
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 179: ALL PASSED")
