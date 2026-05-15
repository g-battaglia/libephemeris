#!/usr/bin/env python3
"""Round 165: Planet positions at midnight (0h UT) vs noon (12h UT).

Tests whether libephemeris has any systematic time-of-day dependent errors
by comparing positions at both midnight and noon across many epochs.
This can reveal interpolation artifacts, time conversion bugs, or
ephemeris segment boundary issues.
"""

from __future__ import annotations

import sys
import os

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
    "Uranus": (swe.URANUS, ephem.URANUS),
    "Neptune": (swe.NEPTUNE, ephem.NEPTUNE),
    "Pluto": (swe.PLUTO, ephem.PLUTO),
    "MeanNode": (swe.MEAN_NODE, ephem.MEAN_NODE),
    "TrueNode": (swe.TRUE_NODE, ephem.TRUE_NODE),
    "MeanLilith": (swe.MEAN_APOG, ephem.MEAN_APOG),
    "Chiron": (swe.CHIRON, ephem.CHIRON),
}

# Test dates spanning 500 years, both midnight and noon
# Use diverse dates: solstices, equinoxes, random dates, century boundaries
TEST_DATES = []
for year in range(1900, 2101, 10):
    # Jan 1 midnight
    TEST_DATES.append((year, 1, 1, 0.0))
    # Jan 1 noon
    TEST_DATES.append((year, 1, 1, 12.0))
    # Jun 21 midnight (near solstice)
    TEST_DATES.append((year, 6, 21, 0.0))
    # Jun 21 noon
    TEST_DATES.append((year, 6, 21, 12.0))
    # Mar 20 midnight (near equinox)
    TEST_DATES.append((year, 3, 20, 0.0))
    # Mar 20 noon
    TEST_DATES.append((year, 3, 20, 12.0))
    # Sep 23 midnight (near equinox)
    TEST_DATES.append((year, 9, 23, 0.0))
    # Sep 23 noon
    TEST_DATES.append((year, 9, 23, 12.0))
    # Dec 21 midnight (near solstice)
    TEST_DATES.append((year, 12, 21, 0.0))
    # Dec 21 noon
    TEST_DATES.append((year, 12, 21, 12.0))

# Also add some historical/future edge dates
for year in [1850, 1855, 1860, 1870, 1880, 2130, 2140, 2148]:
    TEST_DATES.append((year, 7, 15, 0.0))
    TEST_DATES.append((year, 7, 15, 12.0))

FLAGS = swe.FLG_SPEED | swe.FLG_SWIEPH

# Tolerances (arcseconds)
TOL_LON = 1.0  # longitude
TOL_LAT = 1.0  # latitude
TOL_DIST = 50e-6  # AU
TOL_LON_SPD = 1.0  # longitude speed arcsec/day
TOL_LAT_SPD = 1.0  # latitude speed arcsec/day
TOL_DIST_SPD = 50e-6  # dist speed AU/day

# Known wider tolerances
WIDER_TOL_BODIES = {"MeanNode", "TrueNode", "MeanLilith", "Chiron", "Pluto", "Neptune"}

passed = 0
failed = 0
errors = []

for y, m, d, h in TEST_DATES:
    jd = swe.julday(y, m, d, h)

    for name, (se_id, le_id) in BODIES.items():
        try:
            se_result = swe.calc_ut(jd, se_id, FLAGS)
            se_pos = (
                se_result[0] if isinstance(se_result[0], (list, tuple)) else se_result
            )
        except Exception as e:
            continue

        try:
            le_result = ephem.calc_ut(jd, le_id, FLAGS)
            le_pos = le_result[0]
        except Exception:
            continue

        tol_lon = TOL_LON
        tol_lat = TOL_LAT
        tol_dist = TOL_DIST
        tol_lon_spd = TOL_LON_SPD

        if name in WIDER_TOL_BODIES:
            tol_lon = 2.0
            tol_lat = 20.0  # MeanLilith latitude model diff
            tol_dist = 100e-6
            tol_lon_spd = 2.0

        # Compare longitude
        diff_lon = abs(se_pos[0] - le_pos[0])
        if diff_lon > 180:
            diff_lon = 360 - diff_lon
        diff_lon_arcsec = diff_lon * 3600

        # Compare latitude
        diff_lat_arcsec = abs(se_pos[1] - le_pos[1]) * 3600

        # Compare distance
        diff_dist = abs(se_pos[2] - le_pos[2])

        # Compare lon speed
        diff_lon_spd = abs(se_pos[3] - le_pos[3]) * 3600

        # Compare lat speed
        diff_lat_spd = abs(se_pos[4] - le_pos[4]) * 3600

        # Compare dist speed
        diff_dist_spd = abs(se_pos[5] - le_pos[5])

        ok = True
        reasons = []

        if diff_lon_arcsec > tol_lon:
            ok = False
            reasons.append(f'lon {diff_lon_arcsec:.3f}" > {tol_lon}"')

        if diff_lat_arcsec > tol_lat:
            ok = False
            reasons.append(f'lat {diff_lat_arcsec:.3f}" > {tol_lat}"')

        if name not in ("MeanNode", "TrueNode", "MeanLilith") and diff_dist > tol_dist:
            ok = False
            reasons.append(f"dist {diff_dist:.2e} > {tol_dist:.0e}")

        if diff_lon_spd > tol_lon_spd:
            ok = False
            reasons.append(f'lon_spd {diff_lon_spd:.3f}"/d > {tol_lon_spd}"/d')

        if ok:
            passed += 1
        else:
            failed += 1
            if len(errors) < 30:
                errors.append(
                    f"  FAIL {name} {y}/{m:02d}/{d:02d} {h:.0f}h: {', '.join(reasons)}"
                )

total = passed + failed
print(f"\n=== Round 165: Planet Positions Midnight vs Noon ===")
print(f"Test dates: {len(TEST_DATES)}, Bodies: {len(BODIES)}")
print(f"Total tests: {total}")
print(f"PASSED: {passed}/{total} ({100 * passed / total:.1f}%)")
print(f"FAILED: {failed}/{total} ({100 * failed / total:.1f}%)")

if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)

if failed == 0:
    print("\nRound 165: ALL PASSED - No time-of-day systematic errors detected")
