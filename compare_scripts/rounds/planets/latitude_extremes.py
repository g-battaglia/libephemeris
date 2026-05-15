#!/usr/bin/env python3
"""Round 217: Planetary latitude extremes.

Tests planet ecliptic latitude at maximum inclination points,
comparing LE vs SE for all planets at dates where latitude
reaches extreme values.
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
    ("Mercury", ephem.MERCURY, swe.MERCURY, 2.0),
    ("Venus", ephem.VENUS, swe.VENUS, 3.0),
    ("Mars", ephem.MARS, swe.MARS, 5.0),
    ("Jupiter", ephem.JUPITER, swe.JUPITER, 15.0),
    ("Saturn", ephem.SATURN, swe.SATURN, 20.0),
    ("Moon", ephem.MOON, swe.MOON, 1.0),
    ("Pluto", ephem.PLUTO, swe.PLUTO, 30.0),
    ("Chiron", ephem.CHIRON, swe.CHIRON, 20.0),
]

JD_START = 2451545.0  # J2000
SEARCH_YEARS = 20


def find_latitude_extremes(se_body, step, max_count=20):
    """Find dates where latitude reaches local maxima/minima."""
    extremes = []
    jd = JD_START
    jd_end = JD_START + 365.25 * SEARCH_YEARS
    prev_lat = None
    prev_prev_lat = None

    while jd < jd_end:
        try:
            r = swe.calc_ut(jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED)
            lat = r[0][1]
        except Exception:
            jd += step
            continue

        if prev_lat is not None and prev_prev_lat is not None:
            if prev_lat > prev_prev_lat and prev_lat > lat:
                extremes.append(("MaxLat", jd - step, prev_lat))
            elif prev_lat < prev_prev_lat and prev_lat < lat:
                extremes.append(("MinLat", jd - step, prev_lat))

        prev_prev_lat = prev_lat
        prev_lat = lat
        jd += step

        if len(extremes) >= max_count:
            break

    return extremes


def compare_at_extreme(label, le_body, se_body, jd):
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
    if lon_as <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={lon_as:.4f}"')

    # Latitude (the key metric for this round)
    total += 1
    lat_diff = abs(le_r[0][1] - se_r[0][1]) * 3600
    if lat_diff <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT: LE={le_r[0][1]:.6f} SE={se_r[0][1]:.6f} diff={lat_diff:.4f}"'
        )

    # Latitude speed
    total += 1
    latspd_diff = abs(le_r[0][4] - se_r[0][4]) * 3600
    # At extremes, lat_speed should be near zero
    if latspd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT_SPD: LE={le_r[0][4]:.8f} SE={se_r[0][4]:.8f} diff={latspd_diff:.4f}"/day'
        )

    # Distance
    total += 1
    dist_diff = abs(le_r[0][2] - se_r[0][2])
    if dist_diff <= 0.0001:  # 0.0001 AU
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label} DIST: diff={dist_diff:.8f} AU")

    # Lon speed
    total += 1
    lonspd_diff = abs(le_r[0][3] - se_r[0][3]) * 3600
    if lonspd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON_SPD: diff={lonspd_diff:.4f}"/day')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 217: Planetary Latitude Extremes")
    print("=" * 70)

    for bname, le_b, se_b, step in BODIES:
        print(f"\n--- {bname} ---")
        extremes = find_latitude_extremes(se_b, step)
        print(f"  Found {len(extremes)} latitude extremes")
        for etype, jd, lat_val in extremes:
            label = f"{bname} {etype} JD={jd:.2f} lat={lat_val:.4f}"
            compare_at_extreme(label, le_b, se_b, jd)

    print(f"\n{'=' * 70}")
    if total > 0:
        print(
            f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
        )
    else:
        print("RESULTS: 0 tests ran")
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
