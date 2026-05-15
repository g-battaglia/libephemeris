#!/usr/bin/env python3
"""Round 164: House cusp consistency across armc methods.

Compare houses_ex vs houses_armc_ex2 to verify internal consistency.
Both should produce identical cusps when given the same ARMC/obliquity.
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

HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal"),
    ("W", "WholeSgn"),
    ("M", "Morinus"),
    ("B", "Alcabitius"),
    ("T", "Polich"),
]


def le_hsys(ch):
    return ord(ch)


test_configs = []
for year in [2000, 2024, 2025]:
    for lat in [-45, -30, 0, 30, 45, 52, 60]:
        for lon in [0, 90, -73.9]:
            jd = swe.julday(year, 6, 15, 12.0)
            test_configs.append((jd, lat, lon, f"{year} lat={lat} lon={lon}"))

TOL = 0.01  # arcsec — should be exact match for same method

passed = failed = errors = total = 0
failures = []

print(f"Round 164: House Cusp Consistency Across ARMC Methods")
print(f"Testing {len(HOUSE_SYSTEMS)} systems x {len(test_configs)} configs")
print("=" * 90)

for jd, lat, lon, label in test_configs:
    for hsys_ch, hsys_name in HOUSE_SYSTEMS:
        try:
            # Method 1: houses_ex (computes ARMC internally)
            le_r1 = ephem.houses_ex(jd, lat, lon, le_hsys(hsys_ch), FLG_SPEED)
            cusps1 = le_r1[0]
            ascmc1 = le_r1[1]
            armc = ascmc1[2]  # ARMC from houses_ex

            # Get obliquity
            ecl_nut = ephem.calc_ut(jd, -1, 0)
            eps = ecl_nut[0][0]  # true obliquity

            # Method 2: houses_armc_ex2 (using ARMC from method 1)
            le_r2 = ephem.houses_armc_ex2(
                armc, lat, eps, le_hsys(hsys_ch), FLG_SPEED
            )
            cusps2 = le_r2[0]

            # Compare cusps
            n_cusps = min(len(cusps1), len(cusps2), 12)
            for i in range(n_cusps):
                total += 1
                diff = abs(cusps1[i] - cusps2[i]) * 3600.0
                if diff > 180 * 3600:
                    diff = 360 * 3600 - diff
                if diff <= TOL:
                    passed += 1
                else:
                    failed += 1
                    msg = f'  FAIL {label} {hsys_name} cusp[{i + 1}]: ex={cusps1[i]:.6f} armc={cusps2[i]:.6f} diff={diff:.4f}"'
                    failures.append(msg)
                    if len(failures) <= 15:
                        print(msg)

        except Exception as e:
            errors += 1
            if errors <= 3:
                print(f"  ERROR {label} {hsys_name}: {e}")

print(f"\n{'=' * 90}")
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\n{len(failures)} failures")
else:
    print("\nAll tests passed!")
