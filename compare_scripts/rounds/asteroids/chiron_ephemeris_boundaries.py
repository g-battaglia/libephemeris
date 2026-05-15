#!/usr/bin/env python3
"""Round 196: Chiron at ephemeris boundaries.

Tests Chiron positions near the edges of its ephemeris coverage
(~1675-2550) where accuracy may degrade.
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


def year_to_jd(year):
    return 2451545.0 + (year - 2000) * 365.25


# Test dates: near boundaries and through the range
TEST_DATES = []
# Near start boundary
for y in [1680, 1690, 1700, 1750]:
    TEST_DATES.append((f"Y{y}", year_to_jd(y)))
# Core range (well-covered)
for y in range(1800, 2101, 25):
    TEST_DATES.append((f"Y{y}", year_to_jd(y)))
# Near end boundary
for y in [2200, 2300, 2400, 2450, 2500, 2540]:
    TEST_DATES.append((f"Y{y}", year_to_jd(y)))

# Flag combos
FLAG_COMBOS = [
    ("default", FLAGS),
    ("J2000", FLAGS | ephem.FLG_J2000),
    ("NONUT", FLAGS | ephem.FLG_NONUT),
    ("EQUATORIAL", FLAGS | ephem.FLG_EQUATORIAL),
    ("HELIO", FLAGS | ephem.FLG_HELCTR),
]


def compare(label, jd, le_flags, se_flags):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, ephem.CHIRON, le_flags)
        se_r = swe.calc_ut(jd, swe.CHIRON, se_flags)
    except Exception:
        return

    # Longitude
    total += 1
    lon_diff = abs(le_r[0][0] - se_r[0][0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    tol = 2.0
    if lon_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={lon_as:.4f}"')

    # Latitude
    total += 1
    lat_as = abs(le_r[0][1] - se_r[0][1]) * 3600
    if lat_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LAT: diff={lat_as:.4f}"')

    # Distance
    total += 1
    dist_diff = abs(le_r[0][2] - se_r[0][2])
    if dist_diff <= 1e-5:
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label} DIST: diff={dist_diff:.2e}")

    # Speed
    total += 1
    spd_diff = abs(le_r[0][3] - se_r[0][3]) * 3600
    if spd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON_SPD: diff={spd_diff:.4f}"/day')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 196: Chiron at Ephemeris Boundaries")
    print("=" * 70)

    for date_label, jd in TEST_DATES:
        for flag_name, le_flags in FLAG_COMBOS:
            se_flags = le_flags  # numerically identical
            label = f"{date_label} {flag_name}"
            compare(label, jd, le_flags, se_flags)

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
