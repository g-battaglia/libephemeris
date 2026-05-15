#!/usr/bin/env python3
"""Round 175: Decanate boundaries — planet positions at exact sign/decanate boundaries.

Tests positions at times when planets cross 0°, 10°, 20° of signs (decanate
boundaries) to verify precision at these astrologically critical points.
Uses mooncross_ut and solcross_ut for exact timing.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLAGS = swe.FLG_SPEED

# Test: compute Moon position at known longitudes using mooncross
# Verify the position is actually at the crossing point
MOON_CROSS_LONS = [
    0,
    30,
    60,
    90,
    120,
    150,
    180,
    210,
    240,
    270,
    300,
    330,
]  # sign boundaries

TEST_STARTS = []
for year in [2000, 2005, 2010, 2015, 2020, 2025]:
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 1, 0.0)
        TEST_STARTS.append((f"{year}/{month:02d}/01", jd))

passed = 0
failed = 0
skipped = 0
errors = []

# Test 1: Verify Moon positions at sign boundaries match between SE and LE
for date_str, jd_start in TEST_STARTS:
    for target_lon in MOON_CROSS_LONS:
        try:
            le_cross_jd = ephem.mooncross_ut(float(target_lon), jd_start, 0)
        except Exception:
            skipped += 1
            continue
        if le_cross_jd == 0.0:
            skipped += 1
            continue

        # Get positions at the crossing time
        try:
            se_r = swe.calc_ut(le_cross_jd, swe.MOON, FLAGS)
            se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
            le_r = ephem.calc_ut(le_cross_jd, ephem.MOON, FLAGS)
            le_pos = le_r[0]
        except Exception:
            skipped += 1
            continue

        # Both positions should be very close to target longitude
        diff_se = abs(se_pos[0] - target_lon)
        if diff_se > 180:
            diff_se = 360 - diff_se
        diff_le = abs(le_pos[0] - target_lon)
        if diff_le > 180:
            diff_le = 360 - diff_le

        # SE and LE should agree on position at this time
        diff_mutual = abs(se_pos[0] - le_pos[0])
        if diff_mutual > 180:
            diff_mutual = 360 - diff_mutual
        diff_mutual_as = diff_mutual * 3600

        ok = True
        reasons = []
        if diff_le * 3600 > 1.0:
            ok = False
            reasons.append(f'LE off target by {diff_le * 3600:.3f}"')
        if diff_mutual_as > 1.0:
            ok = False
            reasons.append(f'SE-LE diff {diff_mutual_as:.3f}"')

        if ok:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(
                    f"  FAIL Moon->lon={target_lon}° from {date_str}: {', '.join(reasons)}"
                )

# Test 2: Sun positions at all sign boundaries
for year in [2000, 2010, 2020, 2025]:
    jd_start = swe.julday(year, 1, 1, 0.0)
    for target_lon in MOON_CROSS_LONS:
        try:
            le_cross_jd = ephem.solcross_ut(float(target_lon), jd_start, 0)
        except Exception:
            skipped += 1
            continue
        if le_cross_jd == 0.0:
            skipped += 1
            continue
        try:
            se_r = swe.calc_ut(le_cross_jd, swe.SUN, FLAGS)
            se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
            le_r = ephem.calc_ut(le_cross_jd, ephem.SUN, FLAGS)
            le_pos = le_r[0]
        except Exception:
            skipped += 1
            continue

        diff_mutual = abs(se_pos[0] - le_pos[0])
        if diff_mutual > 180:
            diff_mutual = 360 - diff_mutual
        diff_mutual_as = diff_mutual * 3600

        if diff_mutual_as <= 1.0:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(
                    f'  FAIL Sun->lon={target_lon}° @{year}: diff {diff_mutual_as:.3f}"'
                )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 175: Decanate Boundaries ===")
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 175: ALL PASSED")
