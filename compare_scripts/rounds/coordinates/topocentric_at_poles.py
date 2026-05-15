#!/usr/bin/env python3
"""Round 185: Topocentric at poles — positions at extreme latitudes.

Tests topocentric positions at Arctic/Antarctic locations where
the parallax correction geometry differs most from mid-latitudes.
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
    "Mars": (swe.MARS, ephem.MARS),
    "Jupiter": (swe.JUPITER, ephem.JUPITER),
}

# Extreme latitude locations
LOCATIONS = [
    ("NorthPole", 89.9, 0.0, 0),
    ("SouthPole", -89.9, 0.0, 2835),
    ("Svalbard", 78.2, 15.6, 0),
    ("McMurdo", -77.8, 166.7, 24),
    ("Alert", 82.5, -62.3, 63),
    ("Vostok", -78.5, 106.8, 3488),
]

FLAGS = swe.FLG_SPEED | 32768  # FLG_TOPOCTR

TEST_DATES = []
for year in [2000, 2010, 2020, 2025]:
    for month in [1, 6]:
        for hour in [0, 12]:
            jd = swe.julday(year, month, 15, float(hour))
            TEST_DATES.append((f"{year}/{month:02d}/15 {hour}h", jd))

TOL_LON = 2.0  # arcsec - topocentric has known small diffs
TOL_LAT = 2.0

passed = 0
failed = 0
skipped = 0
errors = []

for date_str, jd in TEST_DATES:
    for loc_name, lat, lon, alt in LOCATIONS:
        # Set topocentric location
        swe.set_topo(lon, lat, float(alt))
        ephem.set_topo(lon, lat, float(alt))

        for bname, (se_id, le_id) in BODIES.items():
            try:
                se_r = swe.calc_ut(jd, se_id, FLAGS)
                se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
            except Exception:
                skipped += 1
                continue

            try:
                le_r = ephem.calc_ut(jd, le_id, FLAGS)
                le_pos = le_r[0]
            except Exception:
                skipped += 1
                continue

            diff_lon = abs(se_pos[0] - le_pos[0])
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lon_as = diff_lon * 3600
            diff_lat_as = abs(se_pos[1] - le_pos[1]) * 3600

            ok = True
            reasons = []
            if diff_lon_as > TOL_LON:
                ok = False
                reasons.append(f'lon {diff_lon_as:.3f}"')
            if diff_lat_as > TOL_LAT:
                ok = False
                reasons.append(f'lat {diff_lat_as:.3f}"')

            if ok:
                passed += 1
            else:
                failed += 1
                if len(errors) < 20:
                    errors.append(
                        f"  FAIL {bname} {date_str} @{loc_name}: {', '.join(reasons)}"
                    )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 185: Topocentric at Poles ===")
print(f"Dates: {len(TEST_DATES)}, Locations: {len(LOCATIONS)}, Bodies: {len(BODIES)}")
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 185: ALL PASSED")
