#!/usr/bin/env python3
"""Round 100: Full Chart Stress Test — Random Dates/Locations

Simulates real astrological chart calculations at 500 random date/location
combinations. For each, computes all 10 planets + houses and compares SE vs LE.
"""

from __future__ import annotations

import sys
import os
import time
import random

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2

PLANETS = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
    (11, "TrueNode"),
    (12, "MeanLilith"),
]

LON_TOL = 1.5  # arcsec
LAT_TOL = 1.5
CUSP_TOL = 36.0  # arcsec (0.01°)

random.seed(42)  # Reproducible


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 100: Full Chart Stress Test — 500 Random Charts")
    print("=" * 80)

    flags = FLG_SWIEPH | FLG_SPEED

    for chart_num in range(500):
        # Random date: 1900-2100
        jd = random.uniform(2415020.0, 2488069.5)
        lat = random.uniform(-65, 65)  # Avoid polar for Placidus
        lon = random.uniform(-180, 180)

        chart_ok = True

        # --- Planets ---
        for body_id, body_name in PLANETS:
            total += 1
            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception as e:
                errors += 1
                chart_ok = False
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600
            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            # MeanLilith latitude has known ~19" offset
            lat_tol = 25.0 if body_name == "MeanLilith" else LAT_TOL

            if dlon > LON_TOL or dlat > lat_tol:
                failed += 1
                chart_ok = False
                if len(fail_details) < 15:
                    fail_details.append(
                        f"  FAIL [planet] chart#{chart_num} {body_name:12s} JD={jd:.2f}: "
                        f'dLon={dlon:.3f}" dLat={dlat:.3f}"'
                    )
            else:
                passed += 1

        # --- Houses (Placidus) ---
        total += 1
        try:
            se_cusps, se_ascmc = swe.houses(jd, lat, lon, b"P")
            le_result = ephem.houses(jd, lat, lon, ord("P"))
            le_cusps = le_result[0]
            le_ascmc = le_result[1]

            cusps_ok = True
            for i in range(12):
                diff = abs(se_cusps[i] - le_cusps[i]) * 3600
                if diff > 180 * 3600:
                    diff = 360 * 3600 - diff
                if diff > CUSP_TOL:
                    cusps_ok = False
                    break

            if cusps_ok:
                passed += 1
            else:
                failed += 1
                chart_ok = False
        except Exception as e:
            errors += 1

    pct = 100 * passed / max(1, total)
    print(f"\n  Charts tested: 500")
    print(f"  Total checks: {total}")
    print(f"  Passed: {passed} ({pct:.1f}%)")
    print(f"  Failed: {failed}")
    print(f"  Errors: {errors}")

    if fail_details:
        print("\nSample failures:")
        for d in fail_details:
            print(d)

    print()
    print("=" * 80)
    print(
        f"ROUND 100 FINAL: {passed}/{total} passed ({pct:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
