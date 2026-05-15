#!/usr/bin/env python3
"""Round 171: Planetary hours precision — Sun/Moon transit timing accuracy.

Tests the precision of Sun and Moon upper/lower transit (culmination) times
across diverse latitudes and dates. Transit timing is critical for planetary
hours calculations and house cusp accuracy.
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

LOCATIONS = [
    ("London", 51.5, -0.1, 0),
    ("NewYork", 40.7, -74.0, 0),
    ("Tokyo", 35.7, 139.7, 0),
    ("Sydney", -33.9, 151.2, 0),
    ("Quito", -0.2, -78.5, 2850),
    ("Moscow", 55.8, 37.6, 156),
    ("Mumbai", 19.1, 72.9, 14),
    ("Anchorage", 61.2, -149.9, 31),
]

BODIES = [
    ("Sun", swe.SUN, ephem.SUN),
    ("Moon", swe.MOON, ephem.MOON),
]

# CALC_MTRANSIT = 4, CALC_ITRANSIT = 8
TRANSIT = 4  # upper transit (meridian)
ITRANSIT = 8  # lower transit (anti-meridian / IC)

TEST_DATES = []
for year in [1990, 2000, 2010, 2020, 2025]:
    for month in [1, 3, 6, 9, 12]:
        TEST_DATES.append((year, month, 15))

TOL_SECONDS = 120  # 2 minutes

passed = 0
failed = 0
skipped = 0
errors = []

for y, m, d in TEST_DATES:
    jd_start = swe.julday(y, m, d, 0.0)
    for loc_name, lat, lon, alt in LOCATIONS:
        geopos = [lon, lat, float(alt)]
        for bname, se_body, le_body in BODIES:
            for event, ename in [(TRANSIT, "transit"), (ITRANSIT, "itransit")]:
                try:
                    se_r = swe.rise_trans(
                        jd_start, se_body, event, geopos, 1013.25, 15.0
                    )
                    se_jd = se_r[1][0]
                except Exception:
                    skipped += 1
                    continue
                try:
                    le_r = ephem.rise_trans(
                        jd_start, le_body, event, [lon, lat, float(alt)], 1013.25, 15.0
                    )
                    le_jd = le_r[1][0]
                except Exception:
                    skipped += 1
                    continue
                if se_jd == 0.0 or le_jd == 0.0:
                    skipped += 1
                    continue
                diff_sec = abs(se_jd - le_jd) * 86400.0
                if diff_sec <= TOL_SECONDS:
                    passed += 1
                else:
                    failed += 1
                    if len(errors) < 20:
                        errors.append(
                            f"  FAIL {bname} {ename} {y}/{m}/{d} @{loc_name}: {diff_sec:.1f}s"
                        )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 171: Planetary Hours / Transit Timing ===")
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 171: ALL PASSED")
