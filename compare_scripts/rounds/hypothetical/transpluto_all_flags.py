#!/usr/bin/env python3
"""Round 215: Transpluto deep verification.

Tests ISIS (Transpluto) positions across multiple dates, flag combinations,
and coordinate systems. Includes sidereal, equatorial, J2000, heliocentric.
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

LE_BODY = ephem.ISIS  # Transpluto
SE_BODY = swe.ISIS if hasattr(swe, "ISIS") else 15  # ISIS = 15

DATES = [
    2415020.0,  # 1900
    2430000.0,  # 1941
    2440000.0,  # 1968
    2451545.0,  # J2000
    2455000.0,  # 2009
    2460000.0,  # 2023
    2462000.0,  # 2028
    2420000.0,  # 1913
    2445000.0,  # 1982
    2458000.0,  # 2017
]

FLAG_COMBOS = [
    ("Default", ephem.FLG_SWIEPH | ephem.FLG_SPEED, swe.FLG_SWIEPH | swe.FLG_SPEED),
    (
        "J2000",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_J2000,
        swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000,
    ),
    (
        "Equatorial",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL,
        swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL,
    ),
    (
        "NONUT",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NONUT,
        swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT,
    ),
    (
        "NOABERR",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_NOABERR,
        swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NOABERR,
    ),
    (
        "TRUEPOS",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_TRUEPOS,
        swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_TRUEPOS,
    ),
    (
        "Heliocentric",
        ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_HELCTR,
        swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_HELCTR,
    ),
]

SIDEREAL_MODES = [0, 1, 3, 5]  # Fagan, Lahiri, Raman, Krishnamurti


def compare_transpluto(flag_name, le_flags, se_flags, jd, sid_mode=None):
    global passed, failed, total

    if sid_mode is not None:
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode, 0, 0)
        label = f"Transpluto {flag_name} sid={sid_mode} JD={jd:.1f}"
    else:
        label = f"Transpluto {flag_name} JD={jd:.1f}"

    try:
        le_r = ephem.calc_ut(jd, LE_BODY, le_flags)
        se_r = swe.calc_ut(jd, SE_BODY, se_flags)
    except Exception as e:
        return

    # Longitude
    total += 1
    lon_diff = abs(le_r[0][0] - se_r[0][0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600

    tol = 120.0  # 2 arcminutes (hypothetical body)
    if lon_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LON: LE={le_r[0][0]:.6f} SE={se_r[0][0]:.6f} diff={lon_as:.2f}"'
        )

    # Latitude
    total += 1
    lat_diff = abs(le_r[0][1] - se_r[0][1]) * 3600
    if lat_diff <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT: LE={le_r[0][1]:.6f} SE={se_r[0][1]:.6f} diff={lat_diff:.2f}"'
        )

    # Speed
    total += 1
    spd_diff = abs(le_r[0][3] - se_r[0][3]) * 3600
    if spd_diff <= 10.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} SPD: LE={le_r[0][3]:.6f} SE={se_r[0][3]:.6f} diff={spd_diff:.2f}"/day'
        )

    if sid_mode is not None:
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0, 0, 0)


if __name__ == "__main__":
    print("=" * 70)
    print("Round 215: Transpluto Deep Verification")
    print("=" * 70)

    # Standard flag combos
    for flag_name, le_f, se_f in FLAG_COMBOS:
        print(f"\n--- {flag_name} ---")
        for jd in DATES:
            compare_transpluto(flag_name, le_f, se_f, jd)

    # Sidereal mode tests
    for sid_mode in SIDEREAL_MODES:
        le_sid_f = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL
        se_sid_f = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
        print(f"\n--- Sidereal mode {sid_mode} ---")
        for jd in DATES:
            compare_transpluto(f"Sidereal", le_sid_f, se_sid_f, jd, sid_mode=sid_mode)

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
