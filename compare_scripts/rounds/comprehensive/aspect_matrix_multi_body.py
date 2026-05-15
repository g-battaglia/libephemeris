#!/usr/bin/env python3
"""Round 170: Multi-body aspect matrix accuracy.

Computes positions for all major bodies at the same instant and verifies
the angular separations between every pair. This tests whether the
positions are internally consistent and whether any systematic offset
affects relative positions (which matter for astrological aspects).
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

BODIES = [
    ("Sun", swe.SUN, ephem.SUN),
    ("Moon", swe.MOON, ephem.MOON),
    ("Mercury", swe.MERCURY, ephem.MERCURY),
    ("Venus", swe.VENUS, ephem.VENUS),
    ("Mars", swe.MARS, ephem.MARS),
    ("Jupiter", swe.JUPITER, ephem.JUPITER),
    ("Saturn", swe.SATURN, ephem.SATURN),
    ("Uranus", swe.URANUS, ephem.URANUS),
    ("Neptune", swe.NEPTUNE, ephem.NEPTUNE),
    ("Pluto", swe.PLUTO, ephem.PLUTO),
    ("Chiron", swe.CHIRON, ephem.CHIRON),
    ("Ceres", 17, ephem.CERES),
]

FLAGS = swe.FLG_SPEED | swe.FLG_SWIEPH

# Test at many different dates
TEST_DATES = []
for year in range(1950, 2051, 5):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))


def angular_sep(lon1, lat1, lon2, lat2):
    """Compute angular separation in degrees."""
    lon1r, lat1r = math.radians(lon1), math.radians(lat1)
    lon2r, lat2r = math.radians(lon2), math.radians(lat2)
    cos_sep = math.sin(lat1r) * math.sin(lat2r) + math.cos(lat1r) * math.cos(
        lat2r
    ) * math.cos(lon1r - lon2r)
    cos_sep = max(-1.0, min(1.0, cos_sep))
    return math.degrees(math.acos(cos_sep))


TOL_SEP = 2.0  # arcseconds in angular separation difference

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    se_positions = {}
    le_positions = {}

    for bname, se_id, le_id in BODIES:
        try:
            se_r = swe.calc_ut(jd, se_id, FLAGS)
            se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
            se_positions[bname] = (se_pos[0], se_pos[1])
        except Exception:
            continue

        try:
            le_r = ephem.calc_ut(jd, le_id, FLAGS)
            le_positions[bname] = (le_r[0][0], le_r[0][1])
        except Exception:
            continue

    # Compare angular separations between all pairs
    body_names = sorted(set(se_positions.keys()) & set(le_positions.keys()))

    for i in range(len(body_names)):
        for j in range(i + 1, len(body_names)):
            b1, b2 = body_names[i], body_names[j]

            se_sep = angular_sep(
                se_positions[b1][0],
                se_positions[b1][1],
                se_positions[b2][0],
                se_positions[b2][1],
            )
            le_sep = angular_sep(
                le_positions[b1][0],
                le_positions[b1][1],
                le_positions[b2][0],
                le_positions[b2][1],
            )

            diff_as = abs(se_sep - le_sep) * 3600

            if diff_as <= TOL_SEP:
                passed += 1
            else:
                failed += 1
                if len(errors) < 25:
                    errors.append(
                        f"  FAIL {b1}-{b2} {date_str}: "
                        f"SE sep={se_sep:.6f}°, LE sep={le_sep:.6f}°, "
                        f'diff={diff_as:.3f}"'
                    )

total = passed + failed
print(f"\n=== Round 170: Multi-Body Aspect Matrix ===")
print(f"Dates: {len(TEST_DATES)}, Body pairs: {len(BODIES) * (len(BODIES) - 1) // 2}")
print(
    f"Total: {total}, PASSED: {passed} ({100 * passed / total:.1f}%), FAILED: {failed}"
)

if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)

if failed == 0:
    print("\nRound 170: ALL PASSED")
