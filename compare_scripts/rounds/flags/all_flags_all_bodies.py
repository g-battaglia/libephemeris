#!/usr/bin/env python3
"""Round 190: All flags × all bodies stress test.

Exhaustive test of every major body with every meaningful flag combination.
Tests that no crash occurs and positions match SE within tolerance.
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
errors = 0
failures = []

# All major bodies
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
    ("MeanLilith", ephem.MEAN_APOG, swe.MEAN_APOG),
    ("Chiron", ephem.CHIRON, swe.CHIRON),
    ("Ceres", ephem.CERES, 17),
    ("Pallas", ephem.PALLAS, 18),
]

# Flag combinations to test
FLAG_COMBOS = [
    ("default", 0),
    ("SPEED", ephem.FLG_SPEED),
    ("J2000", ephem.FLG_J2000),
    ("J2000+SPEED", ephem.FLG_J2000 | ephem.FLG_SPEED),
    ("NONUT", ephem.FLG_NONUT),
    ("NONUT+SPEED", ephem.FLG_NONUT | ephem.FLG_SPEED),
    ("EQUATORIAL", ephem.FLG_EQUATORIAL),
    ("EQ+SPEED", ephem.FLG_EQUATORIAL | ephem.FLG_SPEED),
    ("EQ+J2000", ephem.FLG_EQUATORIAL | ephem.FLG_J2000),
    ("EQ+J2000+SPEED", ephem.FLG_EQUATORIAL | ephem.FLG_J2000 | ephem.FLG_SPEED),
    ("NOABERR", ephem.FLG_NOABERR),
    ("NOABERR+SPEED", ephem.FLG_NOABERR | ephem.FLG_SPEED),
    ("TRUEPOS", ephem.FLG_TRUEPOS),
    ("TRUEPOS+SPEED", ephem.FLG_TRUEPOS | ephem.FLG_SPEED),
    ("XYZ", ephem.FLG_XYZ),
    ("XYZ+SPEED", ephem.FLG_XYZ | ephem.FLG_SPEED),
    ("RADIANS", ephem.FLG_RADIANS),
    ("RADIANS+SPEED", ephem.FLG_RADIANS | ephem.FLG_SPEED),
    ("NONUT+NOABERR", ephem.FLG_NONUT | ephem.FLG_NOABERR),
    ("TRUEPOS+EQ", ephem.FLG_TRUEPOS | ephem.FLG_EQUATORIAL),
]

# Skip heliocentric for Moon, nodes, Lilith
HELIO_SKIP = {"Moon", "MeanNode", "TrueNode", "MeanLilith"}

# Test dates
TEST_JDS = [
    2451545.0,  # J2000
    2455197.5,  # 2010 Jan 1
    2458849.5,  # 2020 Jan 1
    2460310.5,  # 2024 Feb 15
]


def compare_pos(
    label, le_body, se_body, jd, le_flags, se_flags, is_xyz=False, is_radians=False
):
    global passed, failed, total, errors

    try:
        le_r = ephem.calc_ut(jd, le_body, le_flags | ephem.FLG_SWIEPH)
    except Exception as e:
        errors += 1
        return

    try:
        se_r = swe.calc_ut(jd, se_body, se_flags | swe.FLG_SWIEPH)
    except Exception as e:
        errors += 1
        return

    # Compare first 3 components
    for idx in range(3):
        total += 1
        le_val = le_r[0][idx]
        se_val = se_r[0][idx]

        if is_xyz:
            diff = abs(le_val - se_val)
            tol = 1e-4  # AU
        elif is_radians and idx < 2:
            diff = abs(le_val - se_val)
            tol = 5e-6  # radians (~1")
        elif idx == 2:
            # Distance
            diff = abs(le_val - se_val)
            tol = 1e-4
        else:
            # Longitude or latitude in degrees
            diff = abs(le_val - se_val)
            if not is_radians and idx == 0 and diff > 180:
                diff = 360 - diff
            diff_arcsec = diff * 3600
            tol_arcsec = 2.0  # generous for stress test
            if diff_arcsec <= tol_arcsec:
                passed += 1
                continue
            else:
                failed += 1
                failures.append(
                    f'  {label} [{idx}]: LE={le_val:.8f} SE={se_val:.8f} diff={diff_arcsec:.4f}"'
                )
                continue

        if diff <= tol:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  {label} [{idx}]: LE={le_val:.8f} SE={se_val:.8f} diff={diff:.2e}"
            )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 190: All Flags × All Bodies Stress Test")
    print("=" * 70)

    for body_name, le_body, se_body in BODIES:
        print(f"\n--- {body_name} ---")

        for flag_name, flag_val in FLAG_COMBOS:
            # Skip heliocentric-like for Moon/nodes
            is_helio = flag_val & ephem.FLG_HELCTR
            if is_helio and body_name in HELIO_SKIP:
                continue

            is_xyz = bool(flag_val & ephem.FLG_XYZ)
            is_radians = bool(flag_val & ephem.FLG_RADIANS)

            le_flags = flag_val
            se_flags = flag_val  # flags are numerically identical

            for jd in TEST_JDS:
                label = f"{body_name} {flag_name} JD={jd:.1f}"
                compare_pos(
                    label, le_body, se_body, jd, le_flags, se_flags, is_xyz, is_radians
                )

    # Also test heliocentric for planets
    print("\n--- Heliocentric ---")
    helio_flag = ephem.FLG_HELCTR | ephem.FLG_SPEED
    for body_name, le_body, se_body in BODIES:
        if body_name in HELIO_SKIP:
            continue
        for jd in TEST_JDS:
            label = f"{body_name} HELIO JD={jd:.1f}"
            compare_pos(
                label, le_body, se_body, jd, helio_flag, helio_flag, False, False
            )

    # Barycentric
    print("--- Barycentric ---")
    bary_flag = ephem.FLG_BARYCTR | ephem.FLG_SPEED
    for body_name, le_body, se_body in BODIES:
        if body_name in {"MeanNode", "TrueNode", "MeanLilith"}:
            continue
        for jd in TEST_JDS:
            label = f"{body_name} BARY JD={jd:.1f}"
            compare_pos(label, le_body, se_body, jd, bary_flag, bary_flag, False, False)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
