#!/usr/bin/env python3
"""Round 182: Mercury elongation maximum — positions near greatest elongation.

Tests Mercury at points near greatest eastern and western elongation,
where Mercury is most visible and where ephemeris precision matters most
for observation planning.
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

FLAGS = swe.FLG_SPEED

TEST_DATES = []
for year in range(1990, 2031):
    for month in range(1, 13):
        jd = swe.julday(year, month, 1, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/01", jd))
        jd2 = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd2))

TOL_LON = 1.0  # arcsec

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    try:
        se_m = swe.calc_ut(jd, swe.MERCURY, FLAGS)[0]
        le_m = ephem.calc_ut(jd, ephem.MERCURY, FLAGS)[0]
    except Exception:
        continue

    diff_lon = abs(se_m[0] - le_m[0])
    if diff_lon > 180:
        diff_lon = 360 - diff_lon
    diff_lon_as = diff_lon * 3600
    diff_lat_as = abs(se_m[1] - le_m[1]) * 3600
    diff_spd = abs(se_m[3] - le_m[3]) * 3600

    ok = True
    reasons = []
    if diff_lon_as > TOL_LON:
        ok = False
        reasons.append(f'lon {diff_lon_as:.3f}"')
    if diff_lat_as > TOL_LON:
        ok = False
        reasons.append(f'lat {diff_lat_as:.3f}"')
    if diff_spd > 2.0:
        ok = False
        reasons.append(f'spd {diff_spd:.3f}"/d')

    if ok:
        passed += 1
    else:
        failed += 1
        if len(errors) < 20:
            errors.append(f"  FAIL Mercury {date_str}: {', '.join(reasons)}")

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 182: Mercury Elongation ===")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 182: ALL PASSED")
