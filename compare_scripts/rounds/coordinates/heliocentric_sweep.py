#!/usr/bin/env python3
"""Round 223: Planet heliocentric positions sweep.

Deep comparison of heliocentric positions for all planets across
multiple dates with FLG_HELCTR flag.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

LE_FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_HELCTR
SE_FLAGS = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_HELCTR

BODIES = [
    ("Mercury", ephem.MERCURY, swe.MERCURY),
    ("Venus", ephem.VENUS, swe.VENUS),
    ("Earth", ephem.EARTH, swe.EARTH),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
    ("Uranus", ephem.URANUS, swe.URANUS),
    ("Neptune", ephem.NEPTUNE, swe.NEPTUNE),
    ("Pluto", ephem.PLUTO, swe.PLUTO),
    ("Chiron", ephem.CHIRON, swe.CHIRON),
    ("Ceres", ephem.CERES, swe.CERES),
]

DATES = [
    2415020.0,
    2420000.0,
    2425000.0,
    2430000.0,
    2435000.0,
    2440000.0,
    2445000.0,
    2451545.0,
    2455000.0,
    2460000.0,
    2462000.0,
    2465000.0,
    2448000.0,
    2453000.0,
    2458000.0,
]


def compare(bname, le_b, se_b, jd):
    global passed, failed, total
    try:
        le_r = ephem.calc_ut(jd, le_b, LE_FLAGS)
        se_r = swe.calc_ut(jd, se_b, SE_FLAGS)
    except:
        return

    label = f"{bname} JD={jd:.1f}"

    # Longitude
    total += 1
    d = abs(le_r[0][0] - se_r[0][0])
    if d > 180:
        d = 360 - d
    das = d * 3600
    tol = 2.0 if bname == "Mercury" else 1.0
    if das <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={das:.4f}"')

    # Latitude
    total += 1
    d = abs(le_r[0][1] - se_r[0][1]) * 3600
    if d <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LAT: diff={d:.4f}"')

    # Distance
    total += 1
    d = abs(le_r[0][2] - se_r[0][2])
    if d <= 0.0001:
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label} DIST: diff={d:.8f} AU")


if __name__ == "__main__":
    print("=" * 70)
    print("Round 223: Planet Heliocentric Positions Sweep")
    print("=" * 70)
    for bname, le_b, se_b in BODIES:
        print(f"\n--- {bname} ---")
        for jd in DATES:
            compare(bname, le_b, se_b, jd)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
