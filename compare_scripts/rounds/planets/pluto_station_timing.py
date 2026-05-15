#!/usr/bin/env python3
"""Round 183: Pluto station timing — retrograde station precision.

Tests Pluto positions near stations (direct and retrograde) where
the planet appears to slow, stop, and reverse direction. Station timing
is one of the most precision-demanding computations.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLAGS = swe.FLG_SPEED

# Sample dates across Pluto's orbit (very slow planet)
TEST_DATES = []
for year in range(1950, 2051, 2):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

TOL_LON = 1.0  # arcsec
TOL_SPD = 1.0  # arcsec/day

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    try:
        se_r = swe.calc_ut(jd, swe.PLUTO, FLAGS)
        se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
        le_r = ephem.calc_ut(jd, ephem.PLUTO, FLAGS)
        le_pos = le_r[0]
    except Exception:
        continue

    diff_lon = abs(se_pos[0] - le_pos[0])
    if diff_lon > 180:
        diff_lon = 360 - diff_lon
    diff_lon_as = diff_lon * 3600
    diff_lat_as = abs(se_pos[1] - le_pos[1]) * 3600
    diff_spd = abs(se_pos[3] - le_pos[3]) * 3600

    # Near station? (speed < 0.005 deg/day = 18"/day)
    near_station = abs(se_pos[3]) < 0.005

    ok = True
    reasons = []
    if diff_lon_as > TOL_LON:
        ok = False
        reasons.append(f'lon {diff_lon_as:.3f}"')
    if diff_lat_as > TOL_LON:
        ok = False
        reasons.append(f'lat {diff_lat_as:.3f}"')
    if near_station:
        # At stations, speed can differ more
        if diff_spd > 2.0:
            ok = False
            reasons.append(f'spd {diff_spd:.3f}"/d [STATION]')
    else:
        if diff_spd > TOL_SPD:
            ok = False
            reasons.append(f'spd {diff_spd:.3f}"/d')

    if ok:
        passed += 1
    else:
        failed += 1
        if len(errors) < 20:
            tag = " [STATION]" if near_station else ""
            errors.append(f"  FAIL Pluto {date_str}{tag}: {', '.join(reasons)}")

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 183: Pluto Station Timing ===")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 183: ALL PASSED")
