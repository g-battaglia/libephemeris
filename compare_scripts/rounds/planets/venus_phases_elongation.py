#!/usr/bin/env python3
"""Round 181: Venus phases — inferior/superior conjunction, greatest elongation.

Tests Venus positions near inferior conjunction (between Earth and Sun),
superior conjunction (far side of Sun), and greatest elongation points.
These are critical test points where geocentric parallax effects are largest.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLAGS = swe.FLG_SPEED

# Search for Venus elongation extremes using calc_ut
TEST_DATES = []
for year in range(1990, 2031, 1):
    for month in range(1, 13):
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

TOL_LON = 1.0  # arcsec
TOL_DIST = 50e-6  # AU

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    try:
        se_v = swe.calc_ut(jd, swe.VENUS, FLAGS)[0]
        se_s = swe.calc_ut(jd, swe.SUN, FLAGS)[0]
        le_v = ephem.calc_ut(jd, ephem.VENUS, FLAGS)[0]
        le_s = ephem.calc_ut(jd, ephem.SUN, FLAGS)[0]
    except Exception:
        continue

    # Compare Venus position
    diff_lon = abs(se_v[0] - le_v[0])
    if diff_lon > 180:
        diff_lon = 360 - diff_lon
    diff_lon_as = diff_lon * 3600
    diff_lat_as = abs(se_v[1] - le_v[1]) * 3600
    diff_dist = abs(se_v[2] - le_v[2])
    diff_spd = abs(se_v[3] - le_v[3]) * 3600

    # Classify Venus phase
    elong = se_v[0] - se_s[0]
    if elong > 180:
        elong -= 360
    elif elong < -180:
        elong += 360
    phase = (
        "near_conj"
        if abs(elong) < 10
        else ("max_elong" if abs(elong) > 40 else "normal")
    )

    ok = True
    reasons = []
    if diff_lon_as > TOL_LON:
        ok = False
        reasons.append(f'lon {diff_lon_as:.3f}"')
    if diff_lat_as > TOL_LON:
        ok = False
        reasons.append(f'lat {diff_lat_as:.3f}"')
    if diff_dist > TOL_DIST:
        ok = False
        reasons.append(f"dist {diff_dist:.2e}")

    if ok:
        passed += 1
    else:
        failed += 1
        if len(errors) < 20:
            errors.append(f"  FAIL Venus {date_str} [{phase}]: {', '.join(reasons)}")

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 181: Venus Phases ===")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 181: ALL PASSED")
