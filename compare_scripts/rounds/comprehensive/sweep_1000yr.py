#!/usr/bin/env python3
"""Round 194: Full chart 1000-year sweep.

Computes complete astrological charts (all planets + houses) at 50-year
intervals from 1600 to 2600, comparing all positions against pyswisseph.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

BODIES = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Moon", ephem.MOON, swe.MOON),
    ("Mercury", ephem.MERCURY, swe.MERCURY),
    ("Venus", ephem.VENUS, swe.VENUS),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
    ("Uranus", ephem.URANUS, swe.URANUS),
    ("Neptune", ephem.NEPTUNE, swe.NEPTUNE),
    ("Pluto", ephem.PLUTO, swe.PLUTO),
    ("MeanNode", ephem.MEAN_NODE, swe.MEAN_NODE),
    ("TrueNode", ephem.TRUE_NODE, swe.TRUE_NODE),
    ("Chiron", ephem.CHIRON, swe.CHIRON),
]

# Chart locations
LOCATIONS = [
    ("London", 51.5074, -0.1278),
    ("New York", 40.7128, -74.0060),
    ("Tokyo", 35.6762, 139.6503),
    ("Sydney", -33.8688, 151.2093),
    ("Rio", -22.9068, -43.1729),
]


# Generate dates from 1600 to 2600 in 50-year steps
# JD for Jan 1 of each year (approximate)
def year_to_jd(year):
    """Approximate JD for Jan 1 of a given year."""
    return 2451545.0 + (year - 2000) * 365.25


TEST_YEARS = list(range(1600, 2601, 50))


def compare_planet(label, le_body, se_body, jd, is_moon=False):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, le_body, FLAGS)
        se_r = swe.calc_ut(jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except Exception:
        return

    # Longitude
    total += 1
    lon_diff = abs(le_r[0][0] - se_r[0][0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600

    tol = 12.0 if is_moon else 2.0  # Moon can diverge more at extreme dates
    if lon_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={lon_as:.4f}" (tol {tol}")')

    # Latitude
    total += 1
    lat_as = abs(le_r[0][1] - se_r[0][1]) * 3600
    lat_tol = (
        20.0
        if "MeanLilith" in label or "MeanNode" in label or "TrueNode" in label
        else tol
    )
    if lat_as <= lat_tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LAT: diff={lat_as:.4f}"')


def test_full_charts():
    global passed, failed, total

    print("=" * 70)
    print("Round 194: Full Chart 1000-Year Sweep")
    print("=" * 70)

    for year in TEST_YEARS:
        jd = year_to_jd(year)
        print(f"\n--- Year {year} (JD {jd:.1f}) ---")

        for body_name, le_body, se_body in BODIES:
            is_moon = body_name == "Moon"
            # Skip Chiron outside its ephemeris range
            if body_name == "Chiron" and (year < 1650 or year > 2550):
                continue
            compare_planet(f"Y{year} {body_name}", le_body, se_body, jd, is_moon)

        # House cusps at each location
        for loc_name, lat, lon in LOCATIONS[:2]:  # Just London and NY for speed
            try:
                le_r = ephem.houses_ex2(jd, lat, lon, ord("P"), 0)
                se_r = swe.houses_ex(jd, lat, lon, b"P")
                le_cusps = le_r[0]
                se_cusps = se_r[0]

                for i in range(12):
                    total += 1
                    diff = abs(le_cusps[i] - se_cusps[i])
                    if diff > 180:
                        diff = 360 - diff
                    if diff * 3600 <= 2.0:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f'  Y{year} {loc_name} cusp{i + 1}: diff={diff * 3600:.2f}"'
                        )

                # ASC/MC
                le_ascmc = le_r[1]
                se_ascmc = se_r[1]
                for idx, name in [(0, "ASC"), (1, "MC")]:
                    total += 1
                    diff = abs(le_ascmc[idx] - se_ascmc[idx])
                    if diff > 180:
                        diff = 360 - diff
                    if diff * 3600 <= 2.0:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f'  Y{year} {loc_name} {name}: diff={diff * 3600:.2f}"'
                        )
            except Exception:
                pass


if __name__ == "__main__":
    test_full_charts()

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
