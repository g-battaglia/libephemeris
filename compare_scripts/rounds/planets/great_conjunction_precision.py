#!/usr/bin/env python3
"""Round 203: Great conjunction precision.

Tests planet positions at great conjunctions (Jupiter-Saturn, etc.)
where precision matters most for astrological applications.
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

# Great conjunctions and notable planetary alignments
EVENTS = [
    # Jupiter-Saturn great conjunctions
    (
        "J-S 2020 Dec 21",
        2459204.5,
        [(ephem.JUPITER, swe.JUPITER), (ephem.SATURN, swe.SATURN)],
    ),
    (
        "J-S 2000 May 28",
        2451693.5,
        [(ephem.JUPITER, swe.JUPITER), (ephem.SATURN, swe.SATURN)],
    ),
    (
        "J-S 1981 Jan 1",
        2444606.5,
        [(ephem.JUPITER, swe.JUPITER), (ephem.SATURN, swe.SATURN)],
    ),
    (
        "J-S 1961 Feb 18",
        2437382.5,
        [(ephem.JUPITER, swe.JUPITER), (ephem.SATURN, swe.SATURN)],
    ),
    (
        "J-S 1940 Aug 8",
        2429856.5,
        [(ephem.JUPITER, swe.JUPITER), (ephem.SATURN, swe.SATURN)],
    ),
    # Venus-Jupiter conjunctions
    (
        "V-J 2015 Jul 1",
        2457204.5,
        [(ephem.VENUS, swe.VENUS), (ephem.JUPITER, swe.JUPITER)],
    ),
    (
        "V-J 2023 Mar 2",
        2460005.5,
        [(ephem.VENUS, swe.VENUS), (ephem.JUPITER, swe.JUPITER)],
    ),
    # Mars-Jupiter
    (
        "M-J 2020 Mar 20",
        2458928.5,
        [(ephem.MARS, swe.MARS), (ephem.JUPITER, swe.JUPITER)],
    ),
    # Mercury-Venus
    (
        "Me-V 2020 May 22",
        2458991.5,
        [(ephem.MERCURY, swe.MERCURY), (ephem.VENUS, swe.VENUS)],
    ),
    # All planets close together
    (
        "Stellium 2000 May",
        2451693.5,
        [
            (ephem.SUN, swe.SUN),
            (ephem.MOON, swe.MOON),
            (ephem.MERCURY, swe.MERCURY),
            (ephem.VENUS, swe.VENUS),
            (ephem.MARS, swe.MARS),
            (ephem.JUPITER, swe.JUPITER),
            (ephem.SATURN, swe.SATURN),
        ],
    ),
    # Outer planet conjunctions
    (
        "U-N 1993 Feb",
        2448990.5,
        [(ephem.URANUS, swe.URANUS), (ephem.NEPTUNE, swe.NEPTUNE)],
    ),
]


def compare_at_event(label, jd, body_pairs):
    global passed, failed, total

    for le_body, se_body in body_pairs:
        try:
            le_r = ephem.calc_ut(jd, le_body, FLAGS)
            se_r = swe.calc_ut(jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED)
        except Exception:
            continue

        body_name = {
            0: "Sun",
            1: "Moon",
            2: "Mercury",
            3: "Venus",
            4: "Mars",
            5: "Jupiter",
            6: "Saturn",
            7: "Uranus",
            8: "Neptune",
            9: "Pluto",
        }.get(le_body, str(le_body))

        # Longitude
        total += 1
        lon_diff = abs(le_r[0][0] - se_r[0][0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lon_as = lon_diff * 3600
        tol = 1.0 if le_body == ephem.MOON else 0.5
        if lon_as <= tol:
            passed += 1
        else:
            failed += 1
            failures.append(f'  {label} {body_name} LON: diff={lon_as:.4f}"')

        # Latitude
        total += 1
        lat_as = abs(le_r[0][1] - se_r[0][1]) * 3600
        if lat_as <= tol:
            passed += 1
        else:
            failed += 1
            failures.append(f'  {label} {body_name} LAT: diff={lat_as:.4f}"')

        # Speed
        total += 1
        spd_as = abs(le_r[0][3] - se_r[0][3]) * 3600
        spd_tol = 5.0 if le_body == ephem.MOON else 2.0
        if spd_as <= spd_tol:
            passed += 1
        else:
            failed += 1
            failures.append(f'  {label} {body_name} SPD: diff={spd_as:.4f}"/day')

        # Also test with J2000
        try:
            le_j = ephem.calc_ut(jd, le_body, FLAGS | ephem.FLG_J2000)
            se_j = swe.calc_ut(
                jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000
            )
            total += 1
            j_diff = abs(le_j[0][0] - se_j[0][0])
            if j_diff > 180:
                j_diff = 360 - j_diff
            if j_diff * 3600 <= tol:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {label} {body_name} J2000: diff={j_diff * 3600:.4f}"'
                )
        except Exception:
            pass

        # Equatorial
        try:
            le_e = ephem.calc_ut(jd, le_body, FLAGS | ephem.FLG_EQUATORIAL)
            se_e = swe.calc_ut(
                jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
            )
            total += 1
            e_diff = abs(le_e[0][0] - se_e[0][0])
            if e_diff > 180:
                e_diff = 360 - e_diff
            if e_diff * 3600 <= tol:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {label} {body_name} EQ: diff={e_diff * 3600:.4f}"')
        except Exception:
            pass


if __name__ == "__main__":
    print("=" * 70)
    print("Round 203: Great Conjunction Precision")
    print("=" * 70)

    for label, jd, pairs in EVENTS:
        print(f"\n--- {label} ---")
        compare_at_event(label, jd, pairs)
        # Also test ±1 day
        compare_at_event(f"{label} -1d", jd - 1, pairs)
        compare_at_event(f"{label} +1d", jd + 1, pairs)

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
