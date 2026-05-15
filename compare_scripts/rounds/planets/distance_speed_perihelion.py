#!/usr/bin/env python3
"""Round 149: Planetary distance speed at perihelion/aphelion.

Compare distance and distance speed (dist_speed) for all planets at dates
near their perihelion/aphelion points, where dist_speed should be near zero
and distance at extremes. Tests precision of orbital mechanics near turning points.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256

BODIES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

# Approximate perihelion/aphelion dates for Earth (affects Sun distance)
# and other planets' closest approach dates (opposition for outer, conjunction for inner)
# We'll sample densely around known perihelion dates

test_cases = []

# Earth perihelion ~Jan 3, aphelion ~Jul 4 each year
for year in range(1950, 2051, 10):
    # Near perihelion
    for day in range(1, 8):
        jd = swe.julday(year, 1, day, 12.0)
        test_cases.append((f"Sun peri {year}-01-{day:02d}", jd, 0))
    # Near aphelion
    for day in range(1, 8):
        jd = swe.julday(year, 7, day, 12.0)
        test_cases.append((f"Sun aph {year}-07-{day:02d}", jd, 0))

# Moon perigee/apogee - cycles every ~27.5 days, sample monthly
for year in [2000, 2010, 2020, 2024, 2025]:
    for month in range(1, 13):
        for day in [1, 8, 15, 22]:
            jd = swe.julday(year, month, day, 12.0)
            test_cases.append((f"Moon {year}-{month:02d}-{day:02d}", jd, 1))

# Inner planets at various phases
for body, bname in [(2, "Mercury"), (3, "Venus")]:
    for year in [1990, 2000, 2010, 2020, 2024]:
        for month in range(1, 13, 2):
            jd = swe.julday(year, month, 15, 12.0)
            test_cases.append((f"{bname} {year}-{month:02d}", jd, body))

# Outer planets - sample less frequently
for body, bname in [
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
]:
    for year in range(1960, 2051, 5):
        jd = swe.julday(year, 6, 15, 12.0)
        test_cases.append((f"{bname} {year}", jd, body))

# Tolerances
TOL_DIST = 1e-6  # AU for distance
TOL_DIST_SPD = 1e-6  # AU/day for distance speed
TOL_LON = 1.0  # arcsec for longitude
TOL_LAT = 1.0  # arcsec for latitude
TOL_LON_SPD = 1.0  # "/day for lon speed

passed = failed = errors = total = 0
failures = []

print(f"Round 149: Planetary Distance Speed at Perihelion/Aphelion")
print(f"Testing {len(test_cases)} cases")
print("=" * 90)

for label, jd, body in test_cases:
    try:
        se_r = swe.calc_ut(jd, body, FLG_SPEED)
        le_r = ephem.calc_ut(jd, body, FLG_SPEED)

        se_data = se_r[0]
        le_data = le_r[0]

        # Compare all 6 components
        comp_names = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]
        tols = [TOL_LON, TOL_LAT, TOL_DIST, TOL_LON_SPD, TOL_LON_SPD, TOL_DIST_SPD]
        units = ['"', '"', "AU", '"/d', '"/d', "AU/d"]

        for i in range(6):
            total += 1
            se_val = se_data[i]
            le_val = le_data[i]

            if i <= 1 or i == 3 or i == 4:  # angular (degrees -> arcsec)
                diff = abs(le_val - se_val) * 3600.0
            else:  # distance (AU)
                diff = abs(le_val - se_val)

            if diff <= tols[i]:
                passed += 1
            else:
                failed += 1
                msg = (
                    f"  FAIL {label} {comp_names[i]}: "
                    f"SE={se_val:.10f} LE={le_val:.10f} diff={diff:.8f}{units[i]}"
                )
                failures.append(msg)
                if len(failures) <= 20:
                    print(msg)
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

print()
print("=" * 90)
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)

if failures:
    print(f"\nTotal failures: {len(failures)}")
    # Categorize by body and component
    cats = {}
    for f in failures:
        for bname in BODIES.values():
            if bname in f:
                for comp in ["dist_spd", "dist", "lon_spd", "lat_spd", "lon", "lat"]:
                    if comp + ":" in f:
                        key = f"{bname} {comp}"
                        cats[key] = cats.get(key, 0) + 1
                        break
    for cat, count in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}")
else:
    print("\nAll tests passed!")
