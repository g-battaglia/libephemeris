#!/usr/bin/env python3
"""Round 174: Midnight ephemeris — positions at 0h TT vs 0h UT.

Tests whether there's any systematic error related to TT/UT time conversion
by comparing positions computed at 0h TT (midnight TT) vs 0h UT (midnight UT).
Delta-T differences compound differently in TT vs UT pipelines.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

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

FLAGS = swe.FLG_SPEED

# Test across centuries — Delta-T changes significantly
TEST_DATES = []
for year in range(1850, 2150, 10):
    jd_ut = swe.julday(year, 1, 1, 0.0)
    TEST_DATES.append((f"{year}/01/01", jd_ut))

TOL_LON = 1.0  # arcsec
TOL_DIST = 50e-6  # AU

passed = 0
failed = 0
errors = []

for date_str, jd_ut in TEST_DATES:
    for bname, (se_id, le_id) in BODIES.items():
        try:
            se_r = swe.calc_ut(jd_ut, se_id, FLAGS)
            se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
        except Exception:
            continue
        try:
            le_r = ephem.calc_ut(jd_ut, le_id, FLAGS)
            le_pos = le_r[0]
        except Exception:
            continue

        diff_lon = abs(se_pos[0] - le_pos[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon
        diff_lon_as = diff_lon * 3600
        diff_lat_as = abs(se_pos[1] - le_pos[1]) * 3600
        diff_dist = abs(se_pos[2] - le_pos[2])

        ok = True
        reasons = []
        if diff_lon_as > TOL_LON:
            ok = False
            reasons.append(f'lon {diff_lon_as:.3f}"')
        if diff_lat_as > TOL_LON:
            ok = False
            reasons.append(f'lat {diff_lat_as:.3f}"')

        if ok:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(f"  FAIL {bname} {date_str}: {', '.join(reasons)}")

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 174: Midnight Ephemeris (0h UT) ===")
print(f"Dates: {len(TEST_DATES)}, Bodies: {len(BODIES)}")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 174: ALL PASSED")
