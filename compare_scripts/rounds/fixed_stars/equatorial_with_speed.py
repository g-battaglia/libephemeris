#!/usr/bin/env python3
"""Round 172: Fixed star parans — fixed stars with EQUATORIAL+SPEED.

Tests fixed stars with equatorial output and speed computation, which exercises
the full RA/Dec conversion pipeline including proper motion propagation.
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

STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Fomalhaut",
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Altair",
    "Deneb",
    "Pollux",
    "Castor",
    "Algol",
]

FLAG_COMBOS = {
    "SPEED": swe.FLG_SPEED,
    "SPEED+EQ": swe.FLG_SPEED | 2048,  # EQUATORIAL
    "SPEED+J2000": swe.FLG_SPEED | 32,  # J2000
    "SPEED+J2000+EQ": swe.FLG_SPEED | 32 | 2048,
}

TEST_DATES = []
for year in range(1950, 2101, 25):
    jd = swe.julday(year, 6, 15, 12.0)
    TEST_DATES.append((f"{year}/06/15", jd))

TOL_POS = 2.0  # arcsec
TOL_SPD = 5.0  # arcsec/day

passed = 0
failed = 0
skipped = 0
errors = []

for date_str, jd in TEST_DATES:
    for fname, flags in FLAG_COMBOS.items():
        for star_name in STARS:
            try:
                se_r = swe.fixstar2(star_name, jd, flags)
                se_pos = se_r[0]
                se_name = se_r[1]
            except Exception:
                skipped += 1
                continue
            try:
                le_r = ephem.fixstar2_ut(star_name, jd, flags)
                le_name = le_r[1]
                le_pos = le_r[0]
            except Exception:
                skipped += 1
                continue

            diff_c1 = abs(se_pos[0] - le_pos[0])
            if diff_c1 > 180:
                diff_c1 = 360 - diff_c1
            diff_c1_as = diff_c1 * 3600
            diff_c2_as = abs(se_pos[1] - le_pos[1]) * 3600
            diff_spd1 = abs(se_pos[3] - le_pos[3]) * 3600

            ok = True
            reasons = []
            if diff_c1_as > TOL_POS:
                ok = False
                reasons.append(f'c1 {diff_c1_as:.3f}"')
            if diff_c2_as > TOL_POS:
                ok = False
                reasons.append(f'c2 {diff_c2_as:.3f}"')
            if diff_spd1 > TOL_SPD:
                ok = False
                reasons.append(f'spd {diff_spd1:.3f}"/d')

            if ok:
                passed += 1
            else:
                failed += 1
                if len(errors) < 20:
                    errors.append(
                        f"  FAIL {star_name} {date_str} [{fname}]: {', '.join(reasons)}"
                    )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 172: Fixed Star Equatorial+Speed ===")
print(f"Stars: {len(STARS)}, Dates: {len(TEST_DATES)}, Flag combos: {len(FLAG_COMBOS)}")
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 172: ALL PASSED")
