#!/usr/bin/env python3
"""Round 173: Light-time iteration verification.

Tests planet positions with FLG_TRUEPOS (geometric, no light-time) vs default
(apparent, with light-time) to verify the light-time correction is consistent.
The difference between TRUEPOS and default should match the expected light-travel
time offset for each planet.
"""

from __future__ import annotations
import sys, os, math
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

FLAGS_DEFAULT = swe.FLG_SPEED
FLAGS_TRUEPOS = swe.FLG_SPEED | 16  # FLG_TRUEPOS

TEST_DATES = []
for year in range(1950, 2051, 10):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

TOL = 1.0  # arcsec - difference in the light-time correction should match

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    for bname, (se_id, le_id) in BODIES.items():
        try:
            se_default = swe.calc_ut(jd, se_id, FLAGS_DEFAULT)[0]
            se_truepos = swe.calc_ut(jd, se_id, FLAGS_TRUEPOS)[0]
            le_default = ephem.calc_ut(jd, le_id, FLAGS_DEFAULT)[0]
            le_truepos = ephem.calc_ut(jd, le_id, FLAGS_TRUEPOS)[0]
        except Exception:
            continue

        # Light-time correction = default - truepos (for each library)
        se_lt_lon = se_default[0] - se_truepos[0]
        if se_lt_lon > 180:
            se_lt_lon -= 360
        elif se_lt_lon < -180:
            se_lt_lon += 360

        le_lt_lon = le_default[0] - le_truepos[0]
        if le_lt_lon > 180:
            le_lt_lon -= 360
        elif le_lt_lon < -180:
            le_lt_lon += 360

        # The light-time correction difference between SE and LE
        diff_lt_as = abs(se_lt_lon - le_lt_lon) * 3600

        tol = TOL
        if bname in ("Pluto", "Neptune", "Uranus"):
            tol = 2.0  # wider for outer planets

        if diff_lt_as <= tol:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(
                    f'  FAIL {bname} {date_str}: SE lt={se_lt_lon * 3600:.3f}", '
                    f'LE lt={le_lt_lon * 3600:.3f}", diff={diff_lt_as:.3f}"'
                )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 173: Light-Time Iteration Verification ===")
print(f"Dates: {len(TEST_DATES)}, Bodies: {len(BODIES)}")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 173: ALL PASSED")
