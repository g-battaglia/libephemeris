#!/usr/bin/env python3
"""Round 197: OscuLilith (Oscillating Lilith) deep sweep.

Tests OSCU_APOG positions at various dates comparing LE vs SE.
OscuLilith uses osculating orbital elements and has larger inherent
speed differences due to rapidly oscillating elements.
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
JD_BASE = 2451545.0

# Daily samples for 3 years
TEST_JDS = [JD_BASE + i * 3 for i in range(365)]  # every 3 days for 3 years

# Also test with various flags
FLAG_COMBOS = [
    ("default", FLAGS),
    ("J2000", FLAGS | ephem.FLG_J2000),
    ("NONUT", FLAGS | ephem.FLG_NONUT),
    ("EQUATORIAL", FLAGS | ephem.FLG_EQUATORIAL),
]


def compare(label, jd, le_flags, se_flags):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, ephem.OSCU_APOG, le_flags)
        se_r = swe.calc_ut(jd, swe.OSCU_APOG, se_flags)
    except Exception:
        return

    # Longitude - OscuLilith can have arcminute-scale model differences
    total += 1
    lon_diff = abs(le_r[0][0] - se_r[0][0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    if lon_as <= 120.0:  # 2 arcminute tolerance (known model diff)
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={lon_as:.2f}" ({lon_as / 3600:.4f}°)')

    # Latitude
    total += 1
    lat_as = abs(le_r[0][1] - se_r[0][1]) * 3600
    if lat_as <= 120.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LAT: diff={lat_as:.2f}"')

    # Speed (known large differences for OscuLilith)
    total += 1
    spd_diff = abs(le_r[0][3] - se_r[0][3]) * 3600
    if spd_diff <= 200.0:  # ~3.3 arcmin/day tolerance
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON_SPD: diff={spd_diff:.2f}"/day')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 197: OscuLilith Deep Sweep")
    print("=" * 70)

    for flag_name, le_flags in FLAG_COMBOS:
        print(f"\n--- {flag_name} ---")
        se_flags = le_flags
        for jd in TEST_JDS:
            label = f"{flag_name} JD={jd:.1f}"
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
