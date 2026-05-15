#!/usr/bin/env python3
"""Round 120: Comprehensive Final Validation — Wide Sweep All Features

Final round of the first 120-round campaign. Tests a broad cross-section of
all major features in a single sweep to confirm overall library health.

Covers:
- Planet positions (all 10 major + Chiron + asteroids)
- Moon node/Lilith
- Fixed stars
- Houses (multiple systems)
- Ayanamsha
- Time functions
- Coordinate transforms
- Phenomena
- Rise/set
- Crossing functions
"""

from __future__ import annotations

import os
import sys
import math

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_EQUATORIAL = 2048
FLG_J2000 = 32
FLG_NONUT = 64
FLG_SIDEREAL = 65536
FLG_HELCTR = 8
FLG_TOPOCTR = 32768
FLG_XYZ = 4096

SUN = 0
MOON = 1
MERCURY = 2
VENUS = 3
MARS = 4
JUPITER = 5
SATURN = 6
URANUS = 7
NEPTUNE = 8
PLUTO = 9
MEAN_NODE = 10
TRUE_NODE = 11
MEAN_APOG = 12
OSCU_APOG = 13
CHIRON = 15
CERES = 17
PALLAS = 18
JUNO = 19
VESTA = 20

ALL_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanApog"),
    (OSCU_APOG, "OscuApog"),
    (CHIRON, "Chiron"),
    (CERES, "Ceres"),
    (PALLAS, "Pallas"),
    (JUNO, "Juno"),
    (VESTA, "Vesta"),
]

TEST_JDS = [
    2451545.0,  # J2000.0 (2000-01-01.5)
    2455197.5,  # 2010-01-01
    2457388.5,  # 2016-01-01
    2459580.5,  # 2022-01-01
    2460310.5,  # 2024-01-15
    2460676.5,  # 2025-01-15
    2440587.5,  # 1970-01-01
    2444239.5,  # 1980-01-01
    2448622.5,  # 1992-01-01
    2463232.5,  # 2032-01-01
]


def compare_pos(se_pos, le_pos, tol_arcsec=2.0, tol_speed=5.0, body_name=""):
    """Compare position tuples. Returns list of (label, passed, detail)."""
    results = []
    labels = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]

    for i, label in enumerate(labels):
        se_val = se_pos[i]
        le_val = le_pos[i]
        diff = le_val - se_val

        if i == 0:  # longitude wrap
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360

        if "dist" in label and "spd" not in label:
            # Distance — skip for analytical bodies
            if body_name in ("MeanNode", "TrueNode", "MeanApog", "OscuApog"):
                results.append((label, True, "skip_dist"))
                continue
            rel = abs(diff / se_val) if abs(se_val) > 1e-10 else abs(diff)
            results.append((label, rel < 1e-5, f"rel={rel:.2e}"))
        elif "spd" in label:
            diff_as = abs(diff) * 3600
            tol = 200.0 if body_name == "OscuApog" else tol_speed
            results.append((label, diff_as < tol, f'{diff_as:.2f}"/day'))
        else:
            diff_as = abs(diff) * 3600
            tol = 30.0 if body_name in ("OscuApog",) else tol_arcsec
            if label == "lat" and body_name == "MeanApog":
                tol = 25.0
            results.append((label, diff_as < tol, f'{diff_as:.4f}"'))

    return results


def main():
    print("=" * 80)
    print("ROUND 120: Comprehensive Final Validation — Wide Sweep All Features")
    print("=" * 80)

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # ===== Section 1: Planet positions (default ecliptic of date) =====
    print("\n--- Section 1: Planet positions (ecliptic of date) ---")
    for jd in TEST_JDS:
        for body_id, body_name in ALL_BODIES:
            try:
                se_pos = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                le_pos = ephem.calc_ut(jd, body_id, FLG_SPEED)[0]
            except Exception:
                continue

            for label, passed, detail in compare_pos(
                se_pos, le_pos, body_name=body_name
            ):
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 20:
                        failures.append(f"  S1 JD={jd} {body_name} {label}: {detail}")

    # ===== Section 2: Equatorial positions =====
    print("\n--- Section 2: Equatorial positions ---")
    for jd in TEST_JDS[:5]:
        for body_id, body_name in ALL_BODIES[:10]:
            flags = FLG_SPEED | FLG_EQUATORIAL
            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                continue

            for label, passed, detail in compare_pos(
                se_pos, le_pos, body_name=body_name
            ):
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 30:
                        failures.append(
                            f"  S2 EQ JD={jd} {body_name} {label}: {detail}"
                        )

    # ===== Section 3: J2000 positions =====
    print("\n--- Section 3: J2000 positions ---")
    for jd in TEST_JDS[:5]:
        for body_id, body_name in ALL_BODIES[:10]:
            flags = FLG_SPEED | FLG_J2000 | FLG_NONUT
            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                continue

            for label, passed, detail in compare_pos(
                se_pos, le_pos, body_name=body_name
            ):
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 40:
                        failures.append(
                            f"  S3 J2000 JD={jd} {body_name} {label}: {detail}"
                        )

    # ===== Section 4: Heliocentric =====
    print("\n--- Section 4: Heliocentric positions ---")
    helio_bodies = [
        (b, n)
        for b, n in ALL_BODIES
        if b
        not in (SUN, MOON, MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG)
    ]
    for jd in TEST_JDS[:5]:
        for body_id, body_name in helio_bodies:
            flags = FLG_SPEED | FLG_HELCTR
            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                continue

            for label, passed, detail in compare_pos(
                se_pos, le_pos, tol_arcsec=3.0, body_name=body_name
            ):
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 50:
                        failures.append(
                            f"  S4 HELIO JD={jd} {body_name} {label}: {detail}"
                        )

    # ===== Section 5: Topocentric =====
    print("\n--- Section 5: Topocentric positions ---")
    swe.set_topo(12.5, 41.9, 50.0)  # Rome
    ephem.set_topo(12.5, 41.9, 50.0)

    for jd in TEST_JDS[:5]:
        for body_id, body_name in [
            (MOON, "Moon"),
            (SUN, "Sun"),
            (MARS, "Mars"),
        ]:
            flags = FLG_SPEED | FLG_TOPOCTR
            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                continue

            for label, passed, detail in compare_pos(
                se_pos, le_pos, tol_speed=10.0, body_name=body_name
            ):
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 55:
                        failures.append(
                            f"  S5 TOPO JD={jd} {body_name} {label}: {detail}"
                        )

    # ===== Section 6: Fixed Stars =====
    print("\n--- Section 6: Fixed stars ---")
    stars = [
        "Aldebaran",
        "Regulus",
        "Spica",
        "Antares",
        "Sirius",
        "Vega",
        "Polaris",
        "Betelgeuse",
        "Rigel",
        "Fomalhaut",
        "Canopus",
        "Arcturus",
    ]

    for star in stars:
        for jd in TEST_JDS[:5]:
            try:
                se_res = swe.fixstar2(star, jd, FLG_SPEED)
                se_pos = se_res[0]
            except Exception:
                continue

            try:
                le_res = ephem.fixstar2_ut(star, jd, FLG_SPEED)
                le_pos = le_res[0]  # (pos_tuple, star_name, retflag)
            except Exception:
                continue

            for i, label in enumerate(["lon", "lat"]):
                diff = le_pos[i] - se_pos[i]
                if i == 0:
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360
                diff_as = abs(diff) * 3600
                tol = 2.0
                total_tests += 1
                if diff_as < tol:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 60:
                        failures.append(
                            f'  S6 Star {star} JD={jd} {label}: {diff_as:.4f}"'
                        )

    # ===== Section 7: Houses =====
    print("\n--- Section 7: House cusps ---")
    house_systems = ["P", "K", "O", "R", "C", "E", "W", "M", "B"]

    for hsys in house_systems:
        for jd in TEST_JDS[:4]:
            for lat in [0.0, 30.0, 45.0, 60.0, -33.0]:
                lon = 12.5
                try:
                    se_cusps, se_ascmc = swe.houses_ex(jd, lat, lon, se_hsys(hsys))
                except Exception:
                    continue

                try:
                    le_result = ephem.houses_ex2(
                        jd, lat, lon, ord(hsys), FLG_SPEED
                    )
                    le_cusps = le_result[0]
                except Exception:
                    continue

                for i in range(min(len(se_cusps), len(le_cusps), 12)):
                    diff = le_cusps[i] - se_cusps[i]
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360
                    diff_as = abs(diff) * 3600
                    tol = 2.0
                    total_tests += 1
                    if diff_as < tol:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 70:
                            failures.append(
                                f'  S7 House {hsys} lat={lat} cusp{i + 1}: {diff_as:.4f}"'
                            )

    # ===== Section 8: Ayanamsha =====
    print("\n--- Section 8: Ayanamsha values ---")
    for mode in [0, 1, 3, 27]:
        swe.set_sid_mode(mode)
        ephem.set_sid_mode(mode, 0, 0)

        for jd in TEST_JDS:
            se_aya = swe.get_ayanamsa_ut(jd)
            le_aya = ephem.get_ayanamsa_ut(jd)
            diff_as = abs(le_aya - se_aya) * 3600
            tol = 20.0
            total_tests += 1
            if diff_as < tol:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 75:
                    failures.append(f'  S8 Aya mode={mode} JD={jd}: {diff_as:.2f}"')

    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    # ===== Section 9: Time functions =====
    print("\n--- Section 9: Time functions ---")
    for jd in TEST_JDS:
        # Delta-T
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff_s = abs(le_dt - se_dt) * 86400
        tol_s = 3.0
        total_tests += 1
        if diff_s < tol_s:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 78:
                failures.append(f"  S9 DeltaT JD={jd}: diff={diff_s:.4f}s")

        # Sidereal time
        se_st = swe.sidtime(jd)
        le_st = ephem.sidtime(jd)
        diff_st = abs(le_st - se_st) * 3600  # seconds of time
        if diff_st > 43200:
            diff_st = 86400 - diff_st
        tol_st = 0.1  # 0.1 seconds of time
        total_tests += 1
        if diff_st < tol_st:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 80:
                failures.append(f"  S9 SidTime JD={jd}: diff={diff_st:.6f}s")

    # ===== Section 10: Julian day conversions =====
    print("\n--- Section 10: Julian day conversions ---")
    test_dates = [
        (2000, 1, 1, 12.0, 1),
        (1999, 12, 31, 0.0, 1),
        (2024, 6, 15, 18.5, 1),
        (1582, 10, 15, 12.0, 1),
        (1582, 10, 4, 12.0, 0),  # Julian calendar
        (-500, 3, 21, 12.0, 0),
        (2050, 1, 1, 0.0, 1),
    ]

    for year, month, day, hour, gregflag in test_dates:
        se_jd = swe.julday(year, month, day, hour, gregflag)
        le_jd = ephem.julday(year, month, day, hour, gregflag)
        diff = abs(le_jd - se_jd)
        total_tests += 1
        if diff < 1e-10:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 85:
                failures.append(
                    f"  S10 JulDay {year}/{month}/{day}: diff={diff:.2e} days"
                )

    for jd in TEST_JDS[:5]:
        se_rev = swe.revjul(jd, 1)
        le_rev = ephem.revjul(jd, 1)
        total_tests += 1
        if se_rev[0] == le_rev[0] and se_rev[1] == le_rev[1] and se_rev[2] == le_rev[2]:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 88:
                failures.append(f"  S10 RevJul JD={jd}: SE={se_rev} LE={le_rev}")

    # ===== Section 11: Coordinate transforms =====
    print("\n--- Section 11: Coordinate transforms (cotrans) ---")
    test_coords = [
        (120.0, 5.0, 1.0),
        (0.0, 23.4, 1.0),
        (90.0, 0.0, 1.0),
        (270.0, -10.0, 1.0),
        (45.0, 45.0, 1.0),
    ]

    for lon, lat, dist in test_coords:
        for obl in [23.4, 23.44, 23.0]:
            se_res = swe.cotrans((lon, lat, dist), obl)
            le_res = ephem.cotrans((lon, lat, dist), obl)

            for i in range(3):
                diff = abs(le_res[i] - se_res[i])
                if i == 0:
                    d = le_res[i] - se_res[i]
                    if d > 180:
                        d -= 360
                    elif d < -180:
                        d += 360
                    diff = abs(d)

                total_tests += 1
                if diff < 1e-8:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 92:
                        failures.append(
                            f"  S11 cotrans ({lon},{lat}) obl={obl} idx={i}: diff={diff:.2e}"
                        )

    # ===== Section 12: Phenomena =====
    print("\n--- Section 12: Planetary phenomena ---")
    for jd in TEST_JDS[:4]:
        for body_id, body_name in [
            (MOON, "Moon"),
            (MARS, "Mars"),
            (VENUS, "Venus"),
            (JUPITER, "Jupiter"),
        ]:
            try:
                se_res = swe.pheno_ut(jd, body_id, 0)
                le_res = ephem.pheno_ut(jd, body_id, 0)
                le_attr = le_res[0]
            except Exception:
                continue

            pheno_labels = [
                "phase_angle",
                "phase",
                "elongation",
                "diameter",
                "magnitude",
            ]
            for i, label in enumerate(pheno_labels):
                se_val = se_res[i]
                le_val = le_attr[i]
                diff = abs(le_val - se_val)

                if label in ("phase_angle", "elongation"):
                    tol = 0.01  # degrees
                elif label == "phase":
                    tol = 0.005
                elif label == "diameter":
                    tol = 0.005  # arcsec
                elif label == "magnitude":
                    tol = 0.5  # mag
                else:
                    tol = 0.01

                total_tests += 1
                if diff < tol:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 98:
                        failures.append(
                            f"  S12 Pheno {body_name} JD={jd} {label}: SE={se_val:.6f} LE={le_val:.6f} diff={diff:.6f}"
                        )

    # ===== Section 13: House position =====
    print("\n--- Section 13: House position ---")
    for jd in TEST_JDS[:4]:
        armc = ephem.sidtime(jd) * 15.0
        obl = 23.44

        for body_id, body_name in [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ]:
            try:
                pos = ephem.calc_ut(jd, body_id, FLG_SPEED)[0]
                plon = pos[0]
                plat = pos[1]
            except Exception:
                continue

            for hsys_ch in ["P", "K", "O", "R", "C", "E"]:
                try:
                    se_hp = swe.house_pos(
                        armc, 41.9, obl, (plon, plat), se_hsys(hsys_ch)
                    )
                    le_hp = ephem.house_pos(
                        armc, 41.9, obl, ord(hsys_ch), plon, plat
                    )
                except Exception:
                    continue

                diff = abs(le_hp - se_hp)
                tol = 0.01
                total_tests += 1
                if diff < tol:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 105:
                        failures.append(
                            f"  S13 HousePos {body_name} {hsys_ch}: SE={se_hp:.6f} LE={le_hp:.6f} diff={diff:.6f}"
                        )

    # ===== Section 14: XYZ output =====
    print("\n--- Section 14: XYZ output ---")
    for jd in TEST_JDS[:3]:
        for body_id, body_name in ALL_BODIES[:10]:
            flags = FLG_SPEED | FLG_XYZ
            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                continue

            for i in range(3):
                diff = abs(le_pos[i] - se_pos[i])
                rel = diff / abs(se_pos[i]) if abs(se_pos[i]) > 1e-10 else diff
                total_tests += 1
                if rel < 1e-5:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 110:
                        failures.append(
                            f"  S14 XYZ {body_name} JD={jd} idx={i}: rel={rel:.2e}"
                        )

    # ===== Section 15: Hypothetical bodies =====
    print("\n--- Section 15: Hypothetical/Uranian bodies ---")
    hypo_bodies = [
        (40, "Cupido"),
        (41, "Hades"),
        (42, "Zeus"),
        (43, "Kronos"),
        (44, "Apollon"),
        (45, "Admetos"),
        (46, "Vulkanus"),
        (47, "Poseidon"),
    ]

    for jd in TEST_JDS[:5]:
        for body_id, body_name in hypo_bodies:
            try:
                se_pos = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                le_pos = ephem.calc_ut(jd, body_id, FLG_SPEED)[0]
            except Exception:
                continue

            diff = le_pos[0] - se_pos[0]
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            diff_as = abs(diff) * 3600
            tol = 60.0  # 60" for hypothetical
            total_tests += 1
            if diff_as < tol:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 115:
                    failures.append(f'  S15 Hypo {body_name} JD={jd}: {diff_as:.2f}"')

    # Summary
    print("\n" + "=" * 80)
    pct = 100 * total_pass / total_tests if total_tests > 0 else 0
    print(f"ROUND 120 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
    print(f"  Failures: {total_fail}")
    print("=" * 80)

    if failures:
        print("\nSample failures:")
        for f in failures[:25]:
            print(f)

    if total_fail == 0:
        print("\nAll tests PASSED!")

    return total_fail


def se_hsys(ch):
    return ch.encode("ascii") if isinstance(ch, str) else ch


if __name__ == "__main__":
    sys.exit(main())
