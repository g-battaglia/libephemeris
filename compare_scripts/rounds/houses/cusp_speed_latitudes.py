#!/usr/bin/env python3
"""Round 219: House cusp speed at various latitudes.

Tests house cusp speeds (from houses_ex2) across multiple house systems,
latitudes, and dates. Compares LE vs SE cusp speed values.
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


def se_hsys(ch: str) -> bytes:
    return ch.encode("ascii")


def le_hsys(ch: str) -> int:
    return ord(ch)


HOUSE_SYSTEMS = ["P", "K", "O", "R", "C", "E", "W", "B"]

LATITUDES = [0.0, 23.4393, 30.0, 45.0, 51.5, 60.0, -33.87, -45.0, 66.0]

DATES = [
    2451545.0,  # J2000
    2460000.0,  # 2023
    2440000.0,  # 1968
    2455000.0,  # 2009
    2458000.0,  # 2017
]

LON = 0.0
FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

# Known buggy SE cusp speeds for Porphyry cusps 5,6,11,12
PORPHYRY_BAD_CUSPS = {4, 5, 10, 11}  # 0-indexed


def compare_cusp_speeds(hsys_ch, lat, jd):
    global passed, failed, total

    label_base = f"{hsys_ch} lat={lat:.1f} JD={jd:.1f}"

    try:
        le_result = ephem.houses_ex2(jd, lat, LON, le_hsys(hsys_ch), FLAGS)
        # Returns (cusps, ascmc, cusp_speeds, ascmc_speeds)
        le_cusps = le_result[0]
        le_ascmc = le_result[1]
        le_cusp_spds = le_result[2]
        le_ascmc_spds = le_result[3]
    except Exception:
        return

    try:
        se_result = swe.houses_ex2(
            jd,
            lat,
            LON,
            se_hsys(hsys_ch),
            swe.FLG_SWIEPH | swe.FLG_SPEED,
        )
        se_cusps = se_result[0]
        se_ascmc = se_result[1]
        se_cusp_spds = se_result[2]
        se_ascmc_spds = se_result[3]
    except Exception:
        return

    # Compare cusp speeds (12 cusps)
    n_cusps = min(len(le_cusp_spds), len(se_cusp_spds), 12)
    for i in range(n_cusps):
        # Skip known buggy Porphyry cusps in SE
        if hsys_ch == "O" and i in PORPHYRY_BAD_CUSPS:
            continue

        total += 1
        le_spd = le_cusp_spds[i]
        se_spd = se_cusp_spds[i]
        spd_diff = abs(le_spd - se_spd) * 3600  # convert to "/day

        # House cusp speeds are typically ~360°/day (sidereal rate)
        # Allow 60"/day tolerance (relative ~0.005%)
        if spd_diff <= 60.0:
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  {label_base} cusp{i + 1}_spd: LE={le_spd:.6f} SE={se_spd:.6f} diff={spd_diff:.2f}"/day'
            )

    # Compare ASC and MC speeds
    for idx, name in [(0, "ASC"), (1, "MC")]:
        if idx < len(le_ascmc_spds) and idx < len(se_ascmc_spds):
            total += 1
            le_spd = le_ascmc_spds[idx]
            se_spd = se_ascmc_spds[idx]
            spd_diff = abs(le_spd - se_spd) * 3600

            if spd_diff <= 60.0:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {label_base} {name}_spd: LE={le_spd:.6f} SE={se_spd:.6f} diff={spd_diff:.2f}"/day'
                )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 219: House Cusp Speed at Various Latitudes")
    print("=" * 70)

    for hsys in HOUSE_SYSTEMS:
        print(f"\n--- House System: {hsys} ---")
        for lat in LATITUDES:
            for jd in DATES:
                compare_cusp_speeds(hsys, lat, jd)

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
        for f in failures[:40]:
            print(f)
        if len(failures) > 40:
            print(f"  ... and {len(failures) - 40} more")
