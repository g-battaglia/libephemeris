#!/usr/bin/env python3
"""Round 202: Alcabitius house_pos deep.

Tests house_pos for Alcabitius and other house systems at various planet
positions and geographic locations.
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
JD = 2451545.0

HOUSE_SYSTEMS = [
    ("Placidus", ord("P"), b"P"),
    ("Koch", ord("K"), b"K"),
    ("Regiomontanus", ord("R"), b"R"),
    ("Campanus", ord("C"), b"C"),
    ("Equal", ord("E"), b"E"),
    ("Porphyry", ord("O"), b"O"),
    ("Alcabitius", ord("B"), b"B"),
]

LOCATIONS = [
    (52.52, 13.405),  # Berlin
    (40.71, -74.01),  # New York
    (35.68, 139.65),  # Tokyo
    (-33.87, 151.21),  # Sydney
    (0.0, 0.0),  # Equator
    (60.17, 24.94),  # Helsinki
    (-45.0, 170.0),  # Southern NZ
]

PLANETS = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Moon", ephem.MOON, swe.MOON),
    ("Mercury", ephem.MERCURY, swe.MERCURY),
    ("Venus", ephem.VENUS, swe.VENUS),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
]


def test_house_pos():
    global passed, failed, total

    print("=" * 70)
    print("Round 202: House Position Deep Test")
    print("=" * 70)

    for lat, lon in LOCATIONS:
        for hsys_name, le_hsys, se_hsys in HOUSE_SYSTEMS:
            # Get ARMC and obliquity
            try:
                le_houses = ephem.houses_ex2(JD, lat, lon, ord("P"), 0)
                armc = le_houses[1][2]  # ARMC
                eps = le_houses[1][4] if len(le_houses[1]) > 4 else 23.4393
            except Exception:
                # Fallback obliquity
                try:
                    ecl = ephem.calc_ut(JD, -1, 0)
                    eps = ecl[0][0]
                except Exception:
                    eps = 23.4393
                continue

            # Get obliquity from ECL_NUT
            try:
                ecl = ephem.calc_ut(JD, -1, 0)
                eps = ecl[0][0]  # true obliquity
            except Exception:
                eps = 23.4393

            for pname, le_body, se_body in PLANETS:
                try:
                    le_pos = ephem.calc_ut(JD, le_body, FLAGS)
                    plon = le_pos[0][0]
                    plat = le_pos[0][1]
                except Exception:
                    continue

                # libephemeris house_pos
                try:
                    le_hp = ephem.house_pos(armc, lat, eps, le_hsys, plon, plat)
                except Exception as e:
                    le_hp = None

                # pyswisseph house_pos
                try:
                    se_hp = swe.house_pos(armc, lat, eps, (plon, plat), se_hsys)
                except Exception:
                    se_hp = None

                if le_hp is None and se_hp is None:
                    continue

                total += 1
                if le_hp is None or se_hp is None:
                    failed += 1
                    failures.append(
                        f"  {hsys_name} lat={lat} {pname}: one returned None"
                    )
                    continue

                diff = abs(le_hp - se_hp)
                if diff > 6:
                    diff = 12 - diff  # wrap around 12 houses

                # Tolerance: 0.01 house units (~0.3° of house)
                if diff <= 0.01:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f"  {hsys_name} lat={lat} {pname}: LE={le_hp:.6f} SE={se_hp:.6f} diff={diff:.6f}"
                    )


if __name__ == "__main__":
    test_house_pos()

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
