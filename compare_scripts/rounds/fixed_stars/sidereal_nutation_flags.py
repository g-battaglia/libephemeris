#!/usr/bin/env python3
"""Round 95: Fixed Stars with Sidereal Mode + Nutation Effects

Tests fixed star positions under:
1. Various sidereal ayanamsha modes (Lahiri, Raman, KP, Fagan-Bradley, etc.)
2. NONUT flag effects on fixed stars
3. J2000 frame fixed stars
4. Combined sidereal + speed flags
5. Nutation amplitude verification on star positions
6. Multiple epochs to check proper motion + sidereal interaction
"""

from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# --- Constants ---
FLG_SPEED = 256
FLG_SIDEREAL = 65536
FLG_NONUT = 64
FLG_J2000 = 32
FLG_EQUATORIAL = 2048
FLG_SWIEPH = 2

SIDM_LAHIRI = 1
SIDM_FAGAN_BRADLEY = 0
SIDM_RAMAN = 3
SIDM_KRISHNAMURTI = 5
SIDM_DELUCE = 2
SIDM_DJWHAL_KHUL = 4
SE_SIDM_YUKTESWAR = 7
SIDM_JN_BHASIN = 8
SIDM_BABYL_KUGLER1 = 9
SIDM_BABYL_KUGLER2 = 10
SIDM_BABYL_KUGLER3 = 11
SIDM_BABYL_HUBER = 12
SIDM_ALDEBARAN_15TAU = 14
SIDM_HIPPARCHOS = 15
SIDM_SASSANIAN = 16
SIDM_GALCENT_0SAG = 17
SIDM_TRUE_CITRA = 27
SIDM_TRUE_REVATI = 28
SIDM_TRUE_PUSHYA = 29
SIDM_TRUE_MULA = 35
SIDM_ARYABHATA = 19

# Stars to test - mix of bright and navigational stars
STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Altair",
    "Pollux",
    "Fomalhaut",
    "Deneb",
    "Castor",
    "Bellatrix",
    "Alnilam",
    "Dubhe",
]

# Ayanamsha modes to test
SIDEREAL_MODES = [
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SIDM_RAMAN, "Raman"),
    (SIDM_KRISHNAMURTI, "KP"),
    (SIDM_DELUCE, "DeLuce"),
    (SE_SIDM_YUKTESWAR, "Yukteswar"),
    (SIDM_JN_BHASIN, "JN_Bhasin"),
    (SIDM_BABYL_KUGLER1, "Babyl_Kugler1"),
    (SIDM_ALDEBARAN_15TAU, "Aldebaran_15Tau"),
    (SIDM_HIPPARCHOS, "Hipparchos"),
    (SIDM_TRUE_CITRA, "True_Citra"),
    (SIDM_TRUE_REVATI, "True_Revati"),
    (SIDM_TRUE_PUSHYA, "True_Pushya"),
]

# Test epochs
EPOCHS = [
    2451545.0,  # J2000.0 (2000-01-01)
    2460000.0,  # 2023
    2440000.0,  # 1968
    2430000.0,  # 1941
    2415020.0,  # 1900-01-01
    2470000.0,  # 2050
]

swe.set_ephe_path(_REF_EPHE_PATH)


def se_fixstar(star, jd, flags):
    """Get fixed star position from pyswisseph."""
    try:
        result = swe.fixstar2(star, jd, flags)
        # result = (pos_tuple, starname, retflag)
        return result[0], result[1]
    except Exception as e:
        return None, str(e)


def le_fixstar(star, jd, flags):
    """Get fixed star position from libephemeris."""
    try:
        result = ephem.fixstar2_ut(star, jd, flags)
        # result = (pos_tuple, starname, retflag)
        return result[0], result[1]
    except Exception as e:
        return None, str(e)


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0

    # Tolerances
    SIDEREAL_TOL = 30.0  # arcsec - known ~14" sidereal offset + star precision
    TROPICAL_TOL = 10.0  # arcsec
    NONUT_TOL = 10.0  # arcsec
    SPEED_TOL = 0.5  # arcsec/day
    LAT_TOL = 5.0  # arcsec for latitude

    print("=" * 80)
    print("ROUND 95: Fixed Stars — Sidereal Mode + Nutation Effects")
    print("=" * 80)

    # =========================================================================
    # PART 1: Sidereal fixed stars across modes and epochs
    # =========================================================================
    print("\n--- PART 1: Sidereal Fixed Stars (multiple modes x stars x epochs) ---")
    p1_pass = 0
    p1_fail = 0
    p1_err = 0

    for sid_mode, sid_name in SIDEREAL_MODES:
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode, 0.0, 0.0)

        for jd in EPOCHS:
            flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL

            for star in STARS:
                total += 1
                se_pos, se_name = se_fixstar(star, jd, flags)
                le_pos, le_name = le_fixstar(star, jd, flags)

                if se_pos is None or le_pos is None:
                    p1_err += 1
                    errors += 1
                    continue

                dlon = abs(se_pos[0] - le_pos[0]) * 3600
                dlat = abs(se_pos[1] - le_pos[1]) * 3600

                # Handle wrap-around
                if dlon > 180 * 3600:
                    dlon = 360 * 3600 - dlon

                if dlon > SIDEREAL_TOL or dlat > LAT_TOL:
                    p1_fail += 1
                    failed += 1
                    if p1_fail <= 10:
                        print(
                            f"  FAIL [{sid_name}] {star} JD={jd:.1f}: "
                            f'dLon={dlon:.2f}" dLat={dlat:.2f}" '
                            f"SE=({se_pos[0]:.6f}, {se_pos[1]:.6f}) "
                            f"LE=({le_pos[0]:.6f}, {le_pos[1]:.6f})"
                        )
                else:
                    p1_pass += 1
                    passed += 1

    print(f"  Part 1 Results: {p1_pass} passed, {p1_fail} failed, {p1_err} errors")

    # =========================================================================
    # PART 2: NONUT flag on fixed stars (tropical)
    # =========================================================================
    print("\n--- PART 2: Fixed Stars with NONUT flag ---")
    p2_pass = 0
    p2_fail = 0
    p2_err = 0

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0.0, 0.0)

    for jd in EPOCHS:
        flags = FLG_SWIEPH | FLG_SPEED | FLG_NONUT

        for star in STARS:
            total += 1
            se_pos, se_name = se_fixstar(star, jd, flags)
            le_pos, le_name = le_fixstar(star, jd, flags)

            if se_pos is None or le_pos is None:
                p2_err += 1
                errors += 1
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            if dlon > NONUT_TOL or dlat > LAT_TOL:
                p2_fail += 1
                failed += 1
                if p2_fail <= 10:
                    print(
                        f"  FAIL [NONUT] {star} JD={jd:.1f}: "
                        f'dLon={dlon:.2f}" dLat={dlat:.2f}" '
                        f"SE=({se_pos[0]:.6f}, {se_pos[1]:.6f}) "
                        f"LE=({le_pos[0]:.6f}, {le_pos[1]:.6f})"
                    )
            else:
                p2_pass += 1
                passed += 1

    print(f"  Part 2 Results: {p2_pass} passed, {p2_fail} failed, {p2_err} errors")

    # =========================================================================
    # PART 3: J2000 frame fixed stars
    # =========================================================================
    print("\n--- PART 3: Fixed Stars in J2000 frame ---")
    p3_pass = 0
    p3_fail = 0
    p3_err = 0

    for jd in EPOCHS:
        flags = FLG_SWIEPH | FLG_SPEED | FLG_J2000 | FLG_NONUT

        for star in STARS:
            total += 1
            se_pos, se_name = se_fixstar(star, jd, flags)
            le_pos, le_name = le_fixstar(star, jd, flags)

            if se_pos is None or le_pos is None:
                p3_err += 1
                errors += 1
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            if dlon > TROPICAL_TOL or dlat > LAT_TOL:
                p3_fail += 1
                failed += 1
                if p3_fail <= 10:
                    print(
                        f"  FAIL [J2000+NONUT] {star} JD={jd:.1f}: "
                        f'dLon={dlon:.2f}" dLat={dlat:.2f}" '
                        f"SE=({se_pos[0]:.6f}, {se_pos[1]:.6f}) "
                        f"LE=({le_pos[0]:.6f}, {le_pos[1]:.6f})"
                    )
            else:
                p3_pass += 1
                passed += 1

    print(f"  Part 3 Results: {p3_pass} passed, {p3_fail} failed, {p3_err} errors")

    # =========================================================================
    # PART 4: Equatorial coordinates (RA/Dec) for fixed stars
    # =========================================================================
    print("\n--- PART 4: Fixed Stars in Equatorial (RA/Dec) ---")
    p4_pass = 0
    p4_fail = 0
    p4_err = 0

    for jd in EPOCHS:
        flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL

        for star in STARS:
            total += 1
            se_pos, se_name = se_fixstar(star, jd, flags)
            le_pos, le_name = le_fixstar(star, jd, flags)

            if se_pos is None or le_pos is None:
                p4_err += 1
                errors += 1
                continue

            # RA in degrees, Dec in degrees
            dra = abs(se_pos[0] - le_pos[0]) * 3600
            ddec = abs(se_pos[1] - le_pos[1]) * 3600

            if dra > 180 * 3600:
                dra = 360 * 3600 - dra

            if dra > TROPICAL_TOL or ddec > LAT_TOL:
                p4_fail += 1
                failed += 1
                if p4_fail <= 10:
                    print(
                        f"  FAIL [EQUATORIAL] {star} JD={jd:.1f}: "
                        f'dRA={dra:.2f}" dDec={ddec:.2f}" '
                        f"SE=({se_pos[0]:.6f}, {se_pos[1]:.6f}) "
                        f"LE=({le_pos[0]:.6f}, {le_pos[1]:.6f})"
                    )
            else:
                p4_pass += 1
                passed += 1

    print(f"  Part 4 Results: {p4_pass} passed, {p4_fail} failed, {p4_err} errors")

    # =========================================================================
    # PART 5: Nutation amplitude consistency check
    # Compute star with and without NONUT, verify the difference matches
    # expected nutation amplitude (~9-17" in longitude)
    # =========================================================================
    print("\n--- PART 5: Nutation Amplitude on Fixed Stars ---")
    p5_pass = 0
    p5_fail = 0
    p5_err = 0

    for jd in EPOCHS:
        for star in STARS[:10]:  # Top 10 stars
            total += 1

            flags_nut = FLG_SWIEPH | FLG_SPEED
            flags_nonut = FLG_SWIEPH | FLG_SPEED | FLG_NONUT

            se_nut, _ = se_fixstar(star, jd, flags_nut)
            se_nonut, _ = se_fixstar(star, jd, flags_nonut)
            le_nut, _ = le_fixstar(star, jd, flags_nut)
            le_nonut, _ = le_fixstar(star, jd, flags_nonut)

            if any(x is None for x in [se_nut, se_nonut, le_nut, le_nonut]):
                p5_err += 1
                errors += 1
                continue

            # Nutation in longitude: difference between with and without NONUT
            se_dpsi = (se_nut[0] - se_nonut[0]) * 3600  # arcsec
            le_dpsi = (le_nut[0] - le_nonut[0]) * 3600

            # The nutation effect should be similar between SE and LE
            dpsi_diff = abs(se_dpsi - le_dpsi)

            # Nutation in obliquity effect on latitude
            se_deps = (se_nut[1] - se_nonut[1]) * 3600
            le_deps = (le_nut[1] - le_nonut[1]) * 3600
            deps_diff = abs(se_deps - le_deps)

            # Tolerance: nutation models should agree within ~1"
            NUT_DIFF_TOL = 2.0  # arcsec

            if dpsi_diff > NUT_DIFF_TOL or deps_diff > NUT_DIFF_TOL:
                p5_fail += 1
                failed += 1
                if p5_fail <= 10:
                    print(
                        f"  FAIL [NUT_AMP] {star} JD={jd:.1f}: "
                        f'SE_dpsi={se_dpsi:.3f}" LE_dpsi={le_dpsi:.3f}" '
                        f'diff={dpsi_diff:.3f}" | '
                        f'SE_deps={se_deps:.3f}" LE_deps={le_deps:.3f}" '
                        f'diff={deps_diff:.3f}"'
                    )
            else:
                p5_pass += 1
                passed += 1

    print(f"  Part 5 Results: {p5_pass} passed, {p5_fail} failed, {p5_err} errors")

    # =========================================================================
    # PART 6: Sidereal + Equatorial combined
    # =========================================================================
    print("\n--- PART 6: Sidereal + Equatorial Combined ---")
    p6_pass = 0
    p6_fail = 0
    p6_err = 0

    MODES_P6 = [
        (SIDM_LAHIRI, "Lahiri"),
        (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
        (SIDM_TRUE_CITRA, "True_Citra"),
    ]

    for sid_mode, sid_name in MODES_P6:
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode, 0.0, 0.0)

        for jd in EPOCHS[:3]:  # J2000, 2023, 1968
            flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL

            for star in STARS[:10]:
                total += 1
                se_pos, _ = se_fixstar(star, jd, flags)
                le_pos, _ = le_fixstar(star, jd, flags)

                if se_pos is None or le_pos is None:
                    p6_err += 1
                    errors += 1
                    continue

                dra = abs(se_pos[0] - le_pos[0]) * 3600
                ddec = abs(se_pos[1] - le_pos[1]) * 3600

                if dra > 180 * 3600:
                    dra = 360 * 3600 - dra

                if dra > SIDEREAL_TOL or ddec > LAT_TOL:
                    p6_fail += 1
                    failed += 1
                    if p6_fail <= 10:
                        print(
                            f"  FAIL [{sid_name}+EQ] {star} JD={jd:.1f}: "
                            f'dRA={dra:.2f}" dDec={ddec:.2f}"'
                        )
                else:
                    p6_pass += 1
                    passed += 1

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0.0, 0.0)

    print(f"  Part 6 Results: {p6_pass} passed, {p6_fail} failed, {p6_err} errors")

    # =========================================================================
    # PART 7: Speed values for fixed stars in sidereal mode
    # =========================================================================
    print("\n--- PART 7: Fixed Star Speeds in Sidereal Mode ---")
    p7_pass = 0
    p7_fail = 0
    p7_err = 0

    MODES_P7 = [
        (SIDM_LAHIRI, "Lahiri"),
        (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    ]

    for sid_mode, sid_name in MODES_P7:
        swe.set_sid_mode(sid_mode)
        ephem.set_sid_mode(sid_mode, 0.0, 0.0)

        for jd in EPOCHS[:3]:
            flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL

            for star in STARS[:10]:
                total += 1
                se_pos, _ = se_fixstar(star, jd, flags)
                le_pos, _ = le_fixstar(star, jd, flags)

                if se_pos is None or le_pos is None:
                    p7_err += 1
                    errors += 1
                    continue

                # Speed in longitude (index 3) and latitude (index 4)
                if len(se_pos) >= 4 and len(le_pos) >= 4:
                    dspd_lon = abs(se_pos[3] - le_pos[3]) * 3600  # arcsec/day
                    dspd_lat = (
                        abs(se_pos[4] - le_pos[4]) * 3600
                        if len(se_pos) >= 5 and len(le_pos) >= 5
                        else 0
                    )

                    if dspd_lon > SPEED_TOL or dspd_lat > SPEED_TOL:
                        p7_fail += 1
                        failed += 1
                        if p7_fail <= 10:
                            print(
                                f"  FAIL [{sid_name}+SPD] {star} JD={jd:.1f}: "
                                f'dSpdLon={dspd_lon:.4f}"/d dSpdLat={dspd_lat:.4f}"/d '
                                f"SE=({se_pos[3]:.8f}) LE=({le_pos[3]:.8f})"
                            )
                    else:
                        p7_pass += 1
                        passed += 1
                else:
                    p7_err += 1
                    errors += 1

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0.0, 0.0)

    print(f"  Part 7 Results: {p7_pass} passed, {p7_fail} failed, {p7_err} errors")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 80)
    print(
        f"ROUND 95 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)

    # Show worst-case analysis
    if failed > 0:
        print('\nNote: Fixed star sidereal mode has a known ~14" systematic offset')
        print(
            "between SE and LE sidereal transformations (documented in accepted differences)."
        )


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
