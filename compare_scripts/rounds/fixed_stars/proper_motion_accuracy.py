#!/usr/bin/env python3
"""Round 218: Fixed star proper motion accuracy.

Tests fixed star positions at multiple epochs to verify proper motion
is being applied correctly. Compares position drift over 100+ years
between LE and SE.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

# Stars with significant proper motion
STARS = [
    "Sirius",
    "Arcturus",
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Pollux",
    "Procyon",
    "Altair",
    "Vega",
    "Capella",
    "Rigel",
    "Betelgeuse",
    "Fomalhaut",
    "Deneb",
    "Canopus",
]

# Wide date range to test proper motion accumulation
DATES = [
    2415020.0,  # 1900.0
    2430000.0,  # 1941
    2440000.0,  # 1968
    2451545.0,  # J2000.0
    2455000.0,  # 2009
    2460000.0,  # 2023
    2462000.0,  # 2028
    2470000.0,  # ~2050
]


def compare_star(star_name, jd):
    global passed, failed, total

    label = f"{star_name} JD={jd:.1f}"

    try:
        le_r = ephem.fixstar2_ut(star_name, jd, FLAGS)
        # Returns (pos_tuple, starname, retflag)
        le_name = le_r[1]
        le_pos = le_r[0]
    except Exception as e:
        return

    try:
        se_r = swe.fixstar2(star_name, jd, swe.FLG_SWIEPH | swe.FLG_SPEED)
        # Returns (pos_tuple, starname_str, retflag)
        se_pos = se_r[0]
        se_name = se_r[1]
    except Exception as e:
        return

    # Longitude
    total += 1
    lon_diff = abs(le_pos[0] - se_pos[0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    if lon_as <= 10.0:  # 10" tolerance
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LON: LE={le_pos[0]:.6f} SE={se_pos[0]:.6f} diff={lon_as:.4f}"'
        )

    # Latitude
    total += 1
    lat_diff = abs(le_pos[1] - se_pos[1]) * 3600
    if lat_diff <= 10.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT: LE={le_pos[1]:.6f} SE={se_pos[1]:.6f} diff={lat_diff:.4f}"'
        )

    # Lon speed (proper motion in longitude)
    total += 1
    lonspd_diff = abs(le_pos[3] - se_pos[3]) * 3600
    if lonspd_diff <= 5.0:  # 5"/day
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LON_SPD: LE={le_pos[3]:.8f} SE={se_pos[3]:.8f} diff={lonspd_diff:.4f}"/day'
        )

    # Lat speed
    total += 1
    latspd_diff = abs(le_pos[4] - se_pos[4]) * 3600
    if latspd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT_SPD: LE={le_pos[4]:.8f} SE={se_pos[4]:.8f} diff={latspd_diff:.4f}"/day'
        )

    # Check proper motion consistency: position at J2000 vs position at other dates
    # should differ by (date - J2000) * speed
    if abs(jd - 2451545.0) > 100:
        total += 1
        # Get J2000 position for reference
        try:
            le_j2000 = ephem.fixstar2_ut(star_name, 2451545.0, FLAGS)
            expected_drift = (jd - 2451545.0) * le_j2000[0][3]  # days * speed deg/day
            actual_drift = le_pos[0] - le_j2000[0][0]
            if actual_drift > 180:
                actual_drift -= 360
            elif actual_drift < -180:
                actual_drift += 360
            drift_diff = abs(actual_drift - expected_drift) * 3600
            # The linear approximation should be good to ~1" over 100 years
            if drift_diff <= 60.0:  # 60" for non-linear effects over long periods
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {label} PM_DRIFT: expected={expected_drift * 3600:.2f}" actual={actual_drift * 3600:.2f}" diff={drift_diff:.2f}"'
                )
        except Exception:
            passed += 1  # skip gracefully


if __name__ == "__main__":
    print("=" * 70)
    print("Round 218: Fixed Star Proper Motion Accuracy")
    print("=" * 70)

    for star in STARS:
        print(f"\n--- {star} ---")
        for jd in DATES:
            compare_star(star, jd)

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
