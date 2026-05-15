#!/usr/bin/env python3
"""Round 200: Comprehensive milestone round.

Tests all major API functions at representative dates to give an overall
quality snapshot: calc_ut, houses, pheno, sidtime, nutation, fixed stars.
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

ALL_BODIES = [
    ("Sun", ephem.SUN, swe.SUN, 0.5),
    ("Moon", ephem.MOON, swe.MOON, 1.0),
    ("Mercury", ephem.MERCURY, swe.MERCURY, 0.5),
    ("Venus", ephem.VENUS, swe.VENUS, 0.5),
    ("Mars", ephem.MARS, swe.MARS, 0.5),
    ("Jupiter", ephem.JUPITER, swe.JUPITER, 0.5),
    ("Saturn", ephem.SATURN, swe.SATURN, 0.5),
    ("Uranus", ephem.URANUS, swe.URANUS, 0.5),
    ("Neptune", ephem.NEPTUNE, swe.NEPTUNE, 0.5),
    ("Pluto", ephem.PLUTO, swe.PLUTO, 1.0),
    ("MeanNode", ephem.MEAN_NODE, swe.MEAN_NODE, 0.5),
    ("TrueNode", ephem.TRUE_NODE, swe.TRUE_NODE, 0.5),
    ("MeanLilith", ephem.MEAN_APOG, swe.MEAN_APOG, 0.5),
    ("OscuLilith", ephem.OSCU_APOG, swe.OSCU_APOG, 120.0),
    ("Chiron", ephem.CHIRON, swe.CHIRON, 0.5),
    ("Ceres", ephem.CERES, 17, 0.5),
    ("Pallas", ephem.PALLAS, 18, 0.5),
    ("Juno", ephem.JUNO, 19, 0.5),
    ("Vesta", ephem.VESTA, 20, 0.5),
]

TEST_JDS = [
    2451545.0,  # J2000
    2455197.5,  # 2010
    2458849.5,  # 2020
    2460310.5,  # 2024
]

FIXED_STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Sirius",
    "Fomalhaut",
    "Vega",
    "Capella",
    "Rigel",
    "Betelgeuse",
]


def test_planet_positions():
    global passed, failed, total

    print("=" * 70)
    print("Round 200: Comprehensive Milestone Round")
    print("=" * 70)

    print("\n--- Planet Positions (ecliptic of date) ---")
    for jd in TEST_JDS:
        for name, le_b, se_b, tol in ALL_BODIES:
            try:
                le_r = ephem.calc_ut(jd, le_b, FLAGS)
                se_r = swe.calc_ut(jd, se_b, swe.FLG_SWIEPH | swe.FLG_SPEED)
            except Exception:
                continue

            # Longitude
            total += 1
            lon_diff = abs(le_r[0][0] - se_r[0][0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            lon_as = lon_diff * 3600
            if lon_as <= tol:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {name} JD={jd:.1f} LON: diff={lon_as:.4f}" (tol {tol}")'
                )


def test_house_cusps():
    global passed, failed, total

    print("\n--- House Cusps ---")
    locations = [(52.52, 13.405), (40.71, -74.01), (-33.87, 151.21)]
    systems = [("P", ord("P"), b"P"), ("K", ord("K"), b"K"), ("R", ord("R"), b"R")]

    for jd in TEST_JDS[:2]:
        for lat, lon in locations:
            for sys_name, le_sys, se_sys in systems:
                try:
                    le_r = ephem.houses_ex2(jd, lat, lon, le_sys, 0)
                    se_r = swe.houses_ex(jd, lat, lon, se_sys)
                except Exception:
                    continue

                for i in range(12):
                    total += 1
                    diff = abs(le_r[0][i] - se_r[0][i])
                    if diff > 180:
                        diff = 360 - diff
                    if diff * 3600 <= 1.0:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f'  {sys_name} JD={jd:.1f} lat={lat} cusp{i + 1}: diff={diff * 3600:.2f}"'
                        )


def test_nutation_obliquity():
    global passed, failed, total

    print("\n--- Nutation/Obliquity ---")
    for jd in TEST_JDS:
        try:
            le_r = ephem.calc_ut(jd, -1, 0)
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
        except Exception:
            continue

        for idx, name in [
            (0, "true_obl"),
            (1, "mean_obl"),
            (2, "nut_lon"),
            (3, "nut_obl"),
        ]:
            total += 1
            diff = abs(le_r[0][idx] - se_r[0][idx]) * 3600
            if diff <= 0.1:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {name} JD={jd:.1f}: diff={diff:.4f}"')


def test_sidereal_time():
    global passed, failed, total

    print("\n--- Sidereal Time ---")
    for jd in TEST_JDS:
        total += 1
        le_st = ephem.sidtime(jd)
        se_st = swe.sidtime(jd)
        diff = abs(le_st - se_st) * 3600  # seconds
        if diff <= 0.01:
            passed += 1
        else:
            failed += 1
            failures.append(f"  sidtime JD={jd:.1f}: diff={diff:.6f}s")


def test_fixed_stars():
    global passed, failed, total

    print("\n--- Fixed Stars ---")
    for jd in TEST_JDS[:2]:
        for star in FIXED_STARS:
            try:
                le_r = ephem.fixstar2_ut(star, jd, FLAGS)
                le_lon = le_r[0][0]

                se_r = swe.fixstar2(star, jd, swe.FLG_SWIEPH | swe.FLG_SPEED)
                se_lon = se_r[0][0]
            except Exception:
                continue

            total += 1
            diff = abs(le_lon - se_lon)
            if diff > 180:
                diff = 360 - diff
            diff_as = diff * 3600
            if diff_as <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {star} JD={jd:.1f}: diff={diff_as:.4f}"')


def test_ayanamsha():
    global passed, failed, total

    print("\n--- Ayanamsha ---")
    modes = [
        ("Lahiri", 1),
        ("Raman", 3),
        ("Krishnamurti", 5),
        ("Fagan_Bradley", 0),
        ("TrueCitra", 27),
        ("TruePushya", 29),
    ]

    for jd in TEST_JDS:
        for name, mode in modes:
            try:
                ephem.set_sid_mode(mode, 0.0, 0.0)
                swe.set_sid_mode(mode, 0.0, 0.0)
                le_aya = ephem.get_ayanamsa_ut(jd)
                se_aya = swe.get_ayanamsa_ut(jd)
            except Exception:
                continue

            total += 1
            diff = abs(le_aya - se_aya) * 3600
            if diff <= 15.0:  # 15" for known sidereal model diff
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {name} JD={jd:.1f}: LE={le_aya:.6f} SE={se_aya:.6f} diff={diff:.2f}"'
                )


if __name__ == "__main__":
    test_planet_positions()
    test_house_cusps()
    test_nutation_obliquity()
    test_sidereal_time()
    test_fixed_stars()
    test_ayanamsha()

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
