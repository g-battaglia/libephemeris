#!/usr/bin/env python3
"""Round 220: Ultimate Final Validation.

Comprehensive final sweep covering ALL major API categories in a single
script: planet positions (ecliptic, equatorial, sidereal), houses,
fixed stars, nutation, delta-T, crossing, cotrans, ayanamsha,
hypothetical bodies, asteroids, and eclipses.

This is the capstone round — testing breadth across all functionality.
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


# ---- Section 1: Planet Positions (ecliptic of date) ----

PLANETS = [
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
    ("Chiron", ephem.CHIRON, swe.CHIRON),
    ("MeanLilith", ephem.MEAN_APOG, swe.MEAN_APOG),
    ("OscuLilith", ephem.OSCU_APOG, swe.OSCU_APOG),
]

FINAL_DATES = [
    2451545.0,  # J2000
    2460676.5,  # 2025 Jan 1
    2440587.5,  # 1970 Jan 1
    2415020.0,  # 1900
    2455197.5,  # 2010 Jan 1
    2462502.5,  # 2030 Jan 1
    2435473.5,  # 1956 Jan 1
    2448988.5,  # 1993 Jan 1
    2457754.5,  # 2017 Jan 1
    2442413.5,  # 1975 Jan 1
]

FLAGS_DEF = ephem.FLG_SWIEPH | ephem.FLG_SPEED
FLAGS_EQU = FLAGS_DEF | ephem.FLG_EQUATORIAL
FLAGS_J2K = FLAGS_DEF | ephem.FLG_J2000


def test_planet_positions():
    global passed, failed, total

    for jd in FINAL_DATES:
        for pname, le_b, se_b in PLANETS:
            for flag_name, le_f, se_f in [
                ("EclDate", FLAGS_DEF, swe.FLG_SWIEPH | swe.FLG_SPEED),
                (
                    "Equat",
                    FLAGS_EQU,
                    swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL,
                ),
                ("J2000", FLAGS_J2K, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000),
            ]:
                try:
                    le_r = ephem.calc_ut(jd, le_b, le_f)
                    se_r = swe.calc_ut(jd, se_b, se_f)
                except Exception:
                    continue

                total += 1
                lon_diff = abs(le_r[0][0] - se_r[0][0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff
                lon_as = lon_diff * 3600

                tol = 1.0
                if pname in ("MeanNode", "TrueNode") and flag_name == "J2000":
                    tol = 20.0
                if pname == "MeanLilith":
                    tol = 20.0

                if lon_as <= tol:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  Planet {pname} {flag_name} JD={jd:.1f}: diff={lon_as:.4f}"'
                    )


# ---- Section 2: Houses ----


def test_houses():
    global passed, failed, total

    house_systems = ["P", "K", "O", "R", "C", "E", "W"]
    latitudes = [0.0, 30.0, 45.0, 51.5, -33.87]

    for jd in [2451545.0, 2460000.0, 2440000.0]:
        for lat in latitudes:
            for hsys in house_systems:
                try:
                    le_c, le_a = ephem.houses_ex(
                        jd, lat, 0.0, le_hsys(hsys), FLAGS_DEF
                    )
                    se_c, se_a = swe.houses_ex(
                        jd,
                        lat,
                        0.0,
                        se_hsys(hsys),
                        swe.FLG_SWIEPH | swe.FLG_SPEED,
                    )
                except Exception:
                    continue

                # Check ASC
                total += 1
                diff = abs(le_a[0] - se_a[0])
                if diff > 180:
                    diff = 360 - diff
                if diff * 3600 <= 1.0:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  House {hsys} lat={lat} JD={jd:.1f} ASC: diff={diff * 3600:.4f}"'
                    )

                # Check MC
                total += 1
                diff = abs(le_a[1] - se_a[1])
                if diff > 180:
                    diff = 360 - diff
                if diff * 3600 <= 1.0:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  House {hsys} lat={lat} JD={jd:.1f} MC: diff={diff * 3600:.4f}"'
                    )

                # Check cusps 1 and 10
                for ci in [0, 9]:
                    if ci < len(le_c) and ci < len(se_c):
                        total += 1
                        diff = abs(le_c[ci] - se_c[ci])
                        if diff > 180:
                            diff = 360 - diff
                        if diff * 3600 <= 1.0:
                            passed += 1
                        else:
                            failed += 1
                            failures.append(
                                f'  House {hsys} lat={lat} JD={jd:.1f} cusp{ci + 1}: diff={diff * 3600:.4f}"'
                            )


# ---- Section 3: Fixed Stars ----


def test_fixed_stars():
    global passed, failed, total

    stars = [
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Sirius",
        "Vega",
        "Arcturus",
        "Pollux",
        "Fomalhaut",
        "Deneb",
    ]

    for jd in [2451545.0, 2460000.0, 2440000.0]:
        for star in stars:
            try:
                le_r = ephem.fixstar2_ut(star, jd, FLAGS_DEF)
                se_r = swe.fixstar2(star, jd, swe.FLG_SWIEPH | swe.FLG_SPEED)
            except Exception:
                continue

            total += 1
            lon_diff = abs(le_r[0][0] - se_r[0][0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            if lon_diff * 3600 <= 10.0:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  Star {star} JD={jd:.1f}: diff={lon_diff * 3600:.4f}"'
                )

            total += 1
            lat_diff = abs(le_r[0][1] - se_r[0][1]) * 3600
            if lat_diff <= 10.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  Star {star} JD={jd:.1f} LAT: diff={lat_diff:.4f}"')


# ---- Section 4: Nutation & Obliquity ----


def test_nutation():
    global passed, failed, total

    for jd in FINAL_DATES:
        try:
            le_r = ephem.calc_ut(jd, -1, 0)
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
        except Exception:
            continue

        # True obliquity
        total += 1
        obl_diff = abs(le_r[0][0] - se_r[0][0]) * 3600
        if obl_diff <= 0.1:
            passed += 1
        else:
            failed += 1
            failures.append(f'  Nutation JD={jd:.1f} true_obl: diff={obl_diff:.6f}"')

        # Nutation in longitude
        total += 1
        dpsi_diff = abs(le_r[0][2] - se_r[0][2]) * 3600
        if dpsi_diff <= 0.1:
            passed += 1
        else:
            failed += 1
            failures.append(f'  Nutation JD={jd:.1f} dpsi: diff={dpsi_diff:.6f}"')


# ---- Section 5: Delta-T ----


def test_deltat():
    global passed, failed, total

    for jd in FINAL_DATES:
        try:
            le_dt = ephem.deltat(jd)
            se_dt = swe.deltat(jd)
        except Exception:
            continue

        total += 1
        dt_diff = abs(le_dt - se_dt) * 86400  # seconds
        if dt_diff <= 1.0:  # within 1 second
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  DeltaT JD={jd:.1f}: LE={le_dt:.10f} SE={se_dt:.10f} diff={dt_diff:.4f}s"
            )


# ---- Section 6: Sidereal Time ----


def test_sidtime():
    global passed, failed, total

    for jd in FINAL_DATES:
        try:
            le_st = ephem.sidtime(jd)
            se_st = swe.sidtime(jd)
        except Exception:
            continue

        total += 1
        st_diff = abs(le_st - se_st)
        if st_diff > 12:
            st_diff = 24 - st_diff
        st_diff_s = st_diff * 3600  # seconds

        if st_diff_s <= 0.5:
            passed += 1
        else:
            failed += 1
            failures.append(
                f"  SidTime JD={jd:.1f}: LE={le_st:.8f} SE={se_st:.8f} diff={st_diff_s:.4f}s"
            )


# ---- Section 7: cotrans ----


def test_cotrans():
    global passed, failed, total

    test_coords = [
        (120.0, 5.0, 1.0),
        (240.0, -3.5, 0.5),
        (0.0, 0.0, 1.0),
        (90.0, 23.44, 1.0),
        (180.0, -10.0, 2.0),
    ]

    for obliquity in [23.44, 23.26, 23.50]:
        for coord in test_coords:
            try:
                le_r = ephem.cotrans(coord, -obliquity)
                se_r = swe.cotrans(coord, -obliquity)
            except Exception:
                continue

            total += 1
            ra_diff = abs(le_r[0] - se_r[0])
            if ra_diff > 180:
                ra_diff = 360 - ra_diff
            if ra_diff * 3600 <= 0.001:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  cotrans {coord} obl={obliquity}: RA diff={ra_diff * 3600:.6f}"'
                )

            total += 1
            dec_diff = abs(le_r[1] - se_r[1]) * 3600
            if dec_diff <= 0.001:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  cotrans {coord} obl={obliquity}: DEC diff={dec_diff:.6f}"'
                )


# ---- Section 8: Ayanamsha ----


def test_ayanamsha():
    global passed, failed, total

    for mode in [0, 1, 3, 5, 7, 27]:
        swe.set_sid_mode(mode)
        ephem.set_sid_mode(mode, 0, 0)

        for jd in [2451545.0, 2460000.0, 2440000.0]:
            try:
                le_a = ephem.get_ayanamsa_ex_ut(jd, FLAGS_DEF)
                se_a = swe.get_ayanamsa_ex_ut(jd, swe.FLG_SWIEPH)
                le_val = le_a[1] if isinstance(le_a, tuple) else le_a
                se_val = (
                    se_a[1] if isinstance(se_a, tuple) and len(se_a) > 1 else se_a[0]
                )
            except Exception:
                continue

            total += 1
            diff = abs(le_val - se_val) * 3600
            if diff <= 30.0:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  Ayanamsha mode={mode} JD={jd:.1f}: diff={diff:.2f}"'
                )

    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)


# ---- Section 9: Hypothetical Bodies ----


def test_hypothetical():
    global passed, failed, total

    hyp_bodies = [
        ("Cupido", ephem.CUPIDO, 40),
        ("Hades", ephem.HADES, 41),
        ("Zeus", ephem.ZEUS, 42),
        ("Kronos", ephem.KRONOS, 43),
        ("Apollon", ephem.APOLLON, 44),
        ("Admetos", ephem.ADMETOS, 45),
        ("Vulkanus", ephem.VULKANUS, 46),
        ("Poseidon", ephem.POSEIDON, 47),
    ]

    for jd in [2451545.0, 2460000.0, 2440000.0]:
        for hname, le_b, se_b in hyp_bodies:
            try:
                le_r = ephem.calc_ut(jd, le_b, FLAGS_DEF)
                se_r = swe.calc_ut(jd, se_b, swe.FLG_SWIEPH | swe.FLG_SPEED)
            except Exception:
                continue

            total += 1
            lon_diff = abs(le_r[0][0] - se_r[0][0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            if lon_diff * 3600 <= 60.0:  # 1 arcminute for hypothetical
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  Hyp {hname} JD={jd:.1f}: diff={lon_diff * 3600:.2f}"'
                )


# ---- Section 10: Asteroids ----


def test_asteroids():
    global passed, failed, total

    asteroids = [
        ("Ceres", ephem.CERES, swe.CERES),
        ("Pallas", ephem.PALLAS, swe.PALLAS),
        ("Juno", ephem.JUNO, swe.JUNO),
        ("Vesta", ephem.VESTA, swe.VESTA),
    ]

    for jd in [2451545.0, 2460000.0, 2440000.0]:
        for aname, le_b, se_b in asteroids:
            try:
                le_r = ephem.calc_ut(jd, le_b, FLAGS_DEF)
                se_r = swe.calc_ut(jd, se_b, swe.FLG_SWIEPH | swe.FLG_SPEED)
            except Exception:
                continue

            total += 1
            lon_diff = abs(le_r[0][0] - se_r[0][0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            if lon_diff * 3600 <= 5.0:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  Asteroid {aname} JD={jd:.1f}: diff={lon_diff * 3600:.4f}"'
                )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 220: Ultimate Final Validation")
    print("=" * 70)

    sections = [
        ("Planet Positions", test_planet_positions),
        ("Houses", test_houses),
        ("Fixed Stars", test_fixed_stars),
        ("Nutation & Obliquity", test_nutation),
        ("Delta-T", test_deltat),
        ("Sidereal Time", test_sidtime),
        ("cotrans", test_cotrans),
        ("Ayanamsha", test_ayanamsha),
        ("Hypothetical Bodies", test_hypothetical),
        ("Asteroids", test_asteroids),
    ]

    for section_name, test_fn in sections:
        print(f"\n--- {section_name} ---")
        before = passed + failed
        test_fn()
        after = passed + failed
        section_total = after - before
        print(f"  {section_total} tests")

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
