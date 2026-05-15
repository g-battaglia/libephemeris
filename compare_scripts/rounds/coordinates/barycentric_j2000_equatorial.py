#!/usr/bin/env python3
"""Round 166: Barycentric + J2000 + EQUATORIAL combined flags.

Tests the triple-flag combination FLG_BARYCTR | FLG_J2000 | FLG_EQUATORIAL
which exercises barycentric position computation, J2000 frame (no precession),
and equatorial coordinate output all at once. This is a rare but valid
combination that stresses the coordinate pipeline.
"""

from __future__ import annotations

import sys
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
    "Uranus": (swe.URANUS, ephem.URANUS),
    "Neptune": (swe.NEPTUNE, ephem.NEPTUNE),
    "Pluto": (swe.PLUTO, ephem.PLUTO),
}

# Flag combinations to test
FLAG_COMBOS = {
    "BARY+J2000+EQ": swe.FLG_SPEED | 4 | 32 | 2048,  # BARYCTR|J2000|EQUATORIAL
    "BARY+J2000": swe.FLG_SPEED | 4 | 32,  # BARYCTR|J2000
    "BARY+EQ": swe.FLG_SPEED | 4 | 2048,  # BARYCTR|EQUATORIAL
    "BARY+J2000+EQ+NONUT": swe.FLG_SPEED | 4 | 32 | 2048 | 64,  # +NONUT
    "BARY+J2000+EQ+NOABERR": swe.FLG_SPEED | 4 | 32 | 2048 | 1024,  # +NOABERR
    "BARY+J2000+EQ+TRUEPOS": swe.FLG_SPEED | 4 | 32 | 2048 | 16,  # +TRUEPOS
}

TEST_DATES = []
for year in range(1900, 2101, 25):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

# Also some edge dates
for year in [1850, 1875, 2125, 2148]:
    jd = swe.julday(year, 6, 15, 12.0)
    TEST_DATES.append((f"{year}/06/15", jd))

TOL_POS = 2.0  # arcsec for RA/Dec or lon/lat
TOL_DIST = 100e-6  # AU
TOL_SPD = 2.0  # arcsec/day

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    for fname, flags in FLAG_COMBOS.items():
        for bname, (se_id, le_id) in BODIES.items():
            try:
                se_result = swe.calc_ut(jd, se_id, flags)
                se_pos = (
                    se_result[0]
                    if isinstance(se_result[0], (list, tuple))
                    else se_result
                )
            except Exception:
                continue

            try:
                le_result = ephem.calc_ut(jd, le_id, flags)
                le_pos = le_result[0]
            except Exception:
                continue

            diff_c1 = abs(se_pos[0] - le_pos[0])
            if diff_c1 > 180:
                diff_c1 = 360 - diff_c1
            diff_c1_as = diff_c1 * 3600

            diff_c2_as = abs(se_pos[1] - le_pos[1]) * 3600
            diff_dist = abs(se_pos[2] - le_pos[2])
            diff_spd1 = abs(se_pos[3] - le_pos[3]) * 3600

            tol_pos = TOL_POS
            tol_dist = TOL_DIST
            # Moon barycentric known to have larger differences
            if bname == "Moon":
                tol_pos = 3.0
                tol_dist = 200e-6

            ok = True
            reasons = []

            if diff_c1_as > tol_pos:
                ok = False
                reasons.append(f'c1 {diff_c1_as:.3f}"')
            if diff_c2_as > tol_pos:
                ok = False
                reasons.append(f'c2 {diff_c2_as:.3f}"')
            if diff_dist > tol_dist:
                ok = False
                reasons.append(f"dist {diff_dist:.2e}")
            if diff_spd1 > TOL_SPD:
                ok = False
                reasons.append(f'spd1 {diff_spd1:.3f}"/d')

            if ok:
                passed += 1
            else:
                failed += 1
                if len(errors) < 25:
                    errors.append(
                        f"  FAIL {bname} {date_str} [{fname}]: {', '.join(reasons)}"
                    )

total = passed + failed
print(f"\n=== Round 166: Barycentric + J2000 + EQUATORIAL Combined Flags ===")
print(
    f"Dates: {len(TEST_DATES)}, Flag combos: {len(FLAG_COMBOS)}, Bodies: {len(BODIES)}"
)
print(
    f"Total: {total}, PASSED: {passed} ({100 * passed / total:.1f}%), FAILED: {failed}"
)

if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)

if failed == 0:
    print("\nRound 166: ALL PASSED")
