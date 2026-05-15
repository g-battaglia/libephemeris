#!/usr/bin/env python3
"""Round 195: Mean node speed perturbations.

Tests the mean lunar node speed computation at fine time intervals,
comparing LE vs SE. LE returns constant analytical mean speed while SE
applies periodic perturbations.
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
JD_BASE = 2451545.0  # J2000


def test_mean_node_speed():
    global passed, failed, total

    print("=" * 70)
    print("Round 195: Mean Node Speed Perturbations")
    print("=" * 70)

    # Test mean node position and speed at daily intervals for 2 years
    print("\n--- Mean Node Daily (2 years) ---")
    max_lon_diff = 0.0
    max_spd_diff = 0.0

    for i in range(0, 730, 1):
        jd = JD_BASE + i

        try:
            le_r = ephem.calc_ut(jd, ephem.MEAN_NODE, FLAGS)
            se_r = swe.calc_ut(jd, swe.MEAN_NODE, swe.FLG_SWIEPH | swe.FLG_SPEED)
        except Exception:
            continue

        # Position
        total += 1
        lon_diff = abs(le_r[0][0] - se_r[0][0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lon_as = lon_diff * 3600
        max_lon_diff = max(max_lon_diff, lon_as)

        if lon_as <= 0.5:
            passed += 1
        else:
            failed += 1
            failures.append(f'  MN day={i} LON: diff={lon_as:.4f}"')

        # Speed
        total += 1
        le_spd = le_r[0][3]
        se_spd = se_r[0][3]
        spd_diff = abs(le_spd - se_spd) * 3600
        max_spd_diff = max(max_spd_diff, spd_diff)

        # Mean node speed: LE returns constant ~-0.0529539°/day
        # SE applies perturbations, so speed varies ~±3"/day
        if spd_diff <= 5.0:  # 5"/day tolerance
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  MN day={i} SPD: LE={le_spd:.8f} SE={se_spd:.8f} diff={spd_diff:.4f}"/day'
            )

    print(f'  Max lon diff: {max_lon_diff:.4f}"')
    print(f'  Max speed diff: {max_spd_diff:.4f}"/day')

    # Test True Node for comparison
    print("\n--- True Node Daily (2 years) ---")
    max_tn_lon = 0.0
    max_tn_spd = 0.0

    for i in range(0, 730, 1):
        jd = JD_BASE + i

        try:
            le_r = ephem.calc_ut(jd, ephem.TRUE_NODE, FLAGS)
            se_r = swe.calc_ut(jd, swe.TRUE_NODE, swe.FLG_SWIEPH | swe.FLG_SPEED)
        except Exception:
            continue

        # Position
        total += 1
        lon_diff = abs(le_r[0][0] - se_r[0][0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lon_as = lon_diff * 3600
        max_tn_lon = max(max_tn_lon, lon_as)

        if lon_as <= 1.0:
            passed += 1
        else:
            failed += 1
            failures.append(f'  TN day={i} LON: diff={lon_as:.4f}"')

        # Speed (True Node has larger speed variation)
        total += 1
        le_spd = le_r[0][3]
        se_spd = se_r[0][3]
        spd_diff = abs(le_spd - se_spd) * 3600
        max_tn_spd = max(max_tn_spd, spd_diff)

        if spd_diff <= 10.0:  # 10"/day for TrueNode speed
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  TN day={i} SPD: LE={le_spd:.8f} SE={se_spd:.8f} diff={spd_diff:.4f}"/day'
            )

    print(f'  Max TN lon diff: {max_tn_lon:.4f}"')
    print(f'  Max TN speed diff: {max_tn_spd:.4f}"/day')

    # Test Mean Lilith
    print("\n--- Mean Lilith Daily (1 year) ---")
    max_ml_lon = 0.0

    for i in range(0, 365, 1):
        jd = JD_BASE + i

        try:
            le_r = ephem.calc_ut(jd, ephem.MEAN_APOG, FLAGS)
            se_r = swe.calc_ut(jd, swe.MEAN_APOG, swe.FLG_SWIEPH | swe.FLG_SPEED)
        except Exception:
            continue

        total += 1
        lon_diff = abs(le_r[0][0] - se_r[0][0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lon_as = lon_diff * 3600
        max_ml_lon = max(max_ml_lon, lon_as)

        if lon_as <= 1.0:
            passed += 1
        else:
            failed += 1
            failures.append(f'  ML day={i} LON: diff={lon_as:.4f}"')

    print(f'  Max ML lon diff: {max_ml_lon:.4f}"')


if __name__ == "__main__":
    test_mean_node_speed()

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
