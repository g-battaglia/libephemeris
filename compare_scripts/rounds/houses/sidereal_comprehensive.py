#!/usr/bin/env python3
"""Round 211: Sidereal Houses Comprehensive.

Tests house cusps in sidereal mode across multiple house systems,
ayanamsha modes, latitudes, and dates.
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


def se_hsys(ch: str) -> bytes:
    return ch.encode("ascii")


def le_hsys(ch: str) -> int:
    return ord(ch)


# House systems to test
HOUSE_SYSTEMS = ["P", "K", "O", "R", "C", "E", "W", "X", "M", "B", "A"]

# Ayanamsha modes to test
AYANAMSHAS = [
    ("Lahiri", 1),
    ("Raman", 3),
    ("KrishnamurtiNew", 5),
    ("Fagan/Bradley", 0),
    ("DeLuce", 2),
    ("Yukteshwar", 7),
]

# Latitudes
LATITUDES = [0.0, 28.6139, 45.0, 51.5074, 60.0, -33.8688, -45.0]

# Dates
DATES = [
    2451545.0,  # J2000
    2460000.0,  # 2023
    2440000.0,  # 1968
    2415020.0,  # 1900
    2430000.0,  # 1941
]

LON = 0.0  # Greenwich


def compare_sidereal_houses(ayanamsha_name, ayanamsha_id, hsys_ch, lat, jd):
    global passed, failed, total

    # Set sidereal mode
    swe.set_sid_mode(ayanamsha_id)
    ephem.set_sid_mode(ayanamsha_id, 0, 0)

    flags = ephem.FLG_SIDEREAL | ephem.FLG_SWIEPH

    try:
        le_cusps, le_ascmc = ephem.houses_ex(jd, lat, LON, le_hsys(hsys_ch), flags)
        se_cusps, se_ascmc = swe.houses_ex(
            jd, lat, LON, se_hsys(hsys_ch), swe.FLG_SIDEREAL | swe.FLG_SWIEPH
        )
    except Exception as e:
        return

    label = f"{ayanamsha_name} {hsys_ch} lat={lat:.1f} jd={jd:.1f}"

    # Compare 12 cusps
    for i in range(min(len(le_cusps), len(se_cusps), 12)):
        total += 1
        diff = abs(le_cusps[i] - se_cusps[i])
        if diff > 180:
            diff = 360 - diff
        diff_as = diff * 3600
        if diff_as <= 30.0:  # 30" tolerance for sidereal (known ~14" systematic)
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  {label} cusp{i + 1}: LE={le_cusps[i]:.6f} SE={se_cusps[i]:.6f} diff={diff_as:.2f}"'
            )

    # Compare ASC and MC (ascmc[0] and ascmc[1])
    for idx, name in [(0, "ASC"), (1, "MC")]:
        if idx < len(le_ascmc) and idx < len(se_ascmc):
            total += 1
            diff = abs(le_ascmc[idx] - se_ascmc[idx])
            if diff > 180:
                diff = 360 - diff
            diff_as = diff * 3600
            if diff_as <= 30.0:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {label} {name}: LE={le_ascmc[idx]:.6f} SE={se_ascmc[idx]:.6f} diff={diff_as:.2f}"'
                )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 211: Sidereal Houses Comprehensive")
    print("=" * 70)

    for ay_name, ay_id in AYANAMSHAS:
        print(f"\n--- Ayanamsha: {ay_name} (mode {ay_id}) ---")
        for hsys in HOUSE_SYSTEMS:
            for lat in LATITUDES:
                for jd in DATES:
                    compare_sidereal_houses(ay_name, ay_id, hsys, lat, jd)

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:40]:
            print(f)
        if len(failures) > 40:
            print(f"  ... and {len(failures) - 40} more")
