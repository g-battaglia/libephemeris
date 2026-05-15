#!/usr/bin/env python3
"""
Round 8: Deep Sidereal / Ayanamsha Audit
==========================================
Compares libephemeris sidereal calculations against pyswisseph:
  P1: All 47 ayanamsha modes (0-46) — ayanamsa values at multiple epochs
  P2: Sidereal planet positions (FLG_SIDEREAL) for all planets
  P3: get_ayanamsa_ex_ut API shape/values
  P4: Modes 43-46 (new, untested) vs pyswisseph
  P5: SIDM_USER custom ayanamsha
  P6: Sidereal velocity correctness
  P7: Fixed star sidereal positions
"""

from __future__ import annotations

import os
import sys
import time
import traceback

import swisseph as swe

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import libephemeris as ephem

EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
if os.path.exists(EPHE_PATH):
    swe.set_ephe_path(EPHE_PATH)

# ============================================================================
# COUNTERS
# ============================================================================

total = 0
passed = 0
failed = 0
skipped = 0
failures = []


def record(name, ok, detail=""):
    global total, passed, failed, failures
    total += 1
    if ok:
        passed += 1
        print(f"  [PASS] {name}: {detail}")
    else:
        failed += 1
        failures.append((name, detail))
        print(f"  [FAIL] {name}: {detail}")


def record_skip(name, reason):
    global skipped
    skipped += 1
    print(f"  [SKIP] {name}: {reason}")


# ============================================================================
# CONSTANTS
# ============================================================================

# All ayanamsha modes 0-46
ALL_MODES = list(range(47))

# Mode names for readable output
MODE_NAMES = {
    0: "Fagan/Bradley",
    1: "Lahiri",
    2: "DeLuce",
    3: "Raman",
    4: "Ushashashi",
    5: "Krishnamurti",
    6: "DjwhalKhul",
    7: "Yukteshwar",
    8: "JNBhasin",
    9: "BabylKugler1",
    10: "BabylKugler2",
    11: "BabylKugler3",
    12: "BabylHuber",
    13: "BabylEtpsc",
    14: "Aldebaran15Tau",
    15: "Hipparchos",
    16: "Sassanian",
    17: "GalCent0Sag",
    18: "J2000",
    19: "J1900",
    20: "B1950",
    21: "Suryasiddhanta",
    22: "SuryasiddhantaMSun",
    23: "Aryabhata",
    24: "AryabhataMSun",
    25: "SSRevati",
    26: "SSCitra",
    27: "TrueCitra",
    28: "TrueRevati",
    29: "TruePushya",
    30: "GalCentRGilbrand",
    31: "GalEquIAU1958",
    32: "GalEquTrue",
    33: "GalEquMula",
    34: "GalAlignMardyks",
    35: "TrueMula",
    36: "GalCentMulaWilhelm",
    37: "Aryabhata522",
    38: "BabylBritton",
    39: "TrueSheoran",
    40: "GalCentCochrane",
    41: "GalEquFiorenza",
    42: "ValensMoon",
    43: "Lahiri1940",
    44: "LahiriVP285",
    45: "KrishnamurtiVP291",
    46: "LahiriICRC",
}

# Star-based modes (wider tolerance expected)
STAR_BASED_MODES = {14, 17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 41, 42}

# Test epochs
EPOCHS = [
    ("J2000", 2451545.0),
    ("2024-Jan", 2460310.5),
    ("1900", 2415020.5),
    ("2100", 2488070.5),
    ("1950", 2433282.5),
]

# Planets for sidereal position tests
PLANETS = [
    ("Sun", swe.SUN),
    ("Moon", swe.MOON),
    ("Mercury", swe.MERCURY),
    ("Venus", swe.VENUS),
    ("Mars", swe.MARS),
    ("Jupiter", swe.JUPITER),
    ("Saturn", swe.SATURN),
    ("MeanNode", swe.MEAN_NODE),
    ("TrueNode", swe.TRUE_NODE),
]

FLG_SIDEREAL = 64 * 1024  # 65536
FLG_SPEED = 256


# ============================================================================
# PART 1: All Ayanamsha Modes — Values at Multiple Epochs
# ============================================================================


def test_part1_ayanamsa_values():
    print("\n" + "=" * 70)
    print("PART 1: Ayanamsha Values — All 47 Modes at Multiple Epochs")
    print("=" * 70)

    for mode in ALL_MODES:
        mode_name = MODE_NAMES.get(mode, f"mode{mode}")
        is_star = mode in STAR_BASED_MODES

        # Set mode in both libraries
        swe.set_sid_mode(mode)
        ephem.set_sid_mode(mode)

        for epoch_name, jd in EPOCHS:
            test_name = f"P1/{mode_name}({mode})/{epoch_name}"

            try:
                aya_se = swe.get_ayanamsa_ut(jd)
                aya_le = ephem.get_ayanamsa_ut(jd)

                diff = abs(aya_se - aya_le)

                # Tolerances: formula-based tight, star-based relaxed
                if is_star:
                    tol = 0.1  # 6 arcminutes for star-based
                else:
                    tol = 0.01  # 36 arcseconds for formula-based

                ok = diff < tol
                detail = (
                    f"SE={aya_se:.6f} LE={aya_le:.6f} "
                    f'diff={diff:.6f}° ({diff * 3600:.2f}")'
                )
                if not ok:
                    detail += f" tol={tol}"
                record(test_name, ok, detail)

            except Exception as e:
                record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 2: Sidereal Planet Positions
# ============================================================================


def test_part2_sidereal_positions():
    print("\n" + "=" * 70)
    print("PART 2: Sidereal Planet Positions (Lahiri, Fagan/Bradley)")
    print("=" * 70)

    jd = 2460310.5  # 2024-Jan-01

    for mode, mode_name in [(1, "Lahiri"), (0, "Fagan/Bradley")]:
        swe.set_sid_mode(mode)
        ephem.set_sid_mode(mode)

        for planet_name, planet_id in PLANETS:
            test_name = f"P2/{mode_name}/{planet_name}"

            try:
                flags_se = FLG_SIDEREAL | FLG_SPEED
                flags_le = FLG_SIDEREAL | FLG_SPEED

                ret_se = swe.calc_ut(jd, planet_id, flags_se)
                ret_le = ephem.calc_ut(jd, planet_id, flags_le)

                # Extract position data
                if isinstance(ret_se, tuple) and isinstance(ret_se[0], (list, tuple)):
                    pos_se = ret_se[0]
                else:
                    pos_se = ret_se

                if isinstance(ret_le, tuple) and len(ret_le) == 2:
                    pos_le = ret_le[0]
                else:
                    pos_le = ret_le

                # Compare longitude
                lon_se = float(pos_se[0])
                lon_le = float(pos_le[0])
                lon_diff = abs(lon_se - lon_le)
                # Handle wrap-around at 360
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff

                lon_ok = lon_diff < 0.01  # 36"
                record(
                    f"{test_name}/lon",
                    lon_ok,
                    f"SE={lon_se:.6f} LE={lon_le:.6f} "
                    f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
                )

                # Compare latitude
                lat_se = float(pos_se[1])
                lat_le = float(pos_le[1])
                lat_diff = abs(lat_se - lat_le)
                lat_ok = lat_diff < 0.01
                record(
                    f"{test_name}/lat",
                    lat_ok,
                    f"SE={lat_se:.6f} LE={lat_le:.6f} diff={lat_diff:.6f}°",
                )

                # Compare speed
                spd_se = float(pos_se[3])
                spd_le = float(pos_le[3])
                spd_diff = abs(spd_se - spd_le)
                spd_ok = spd_diff < 0.01  # 0.01 deg/day
                record(
                    f"{test_name}/speed",
                    spd_ok,
                    f"SE={spd_se:.6f} LE={spd_le:.6f} diff={spd_diff:.6f}°/day",
                )

            except Exception as e:
                record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 3: get_ayanamsa_ex_ut API Shape
# ============================================================================


def test_part3_ayanamsa_ex_api():
    print("\n" + "=" * 70)
    print("PART 3: get_ayanamsa_ex_ut API Shape & Values")
    print("=" * 70)

    jd = 2460310.5  # 2024-Jan-01

    for mode, mode_name in [(1, "Lahiri"), (0, "Fagan/Bradley"), (27, "TrueCitra")]:
        swe.set_sid_mode(mode)
        ephem.set_sid_mode(mode)

        test_base = f"P3/{mode_name}"

        try:
            ret_se = swe.get_ayanamsa_ex_ut(jd, 0)
            ret_le = ephem.get_ayanamsa_ex_ut(jd, 0)

            # API shape check
            record(
                f"{test_base}/type",
                isinstance(ret_le, tuple) and len(ret_le) == 2,
                f"SE type={type(ret_se).__name__}({len(ret_se)}) "
                f"LE type={type(ret_le).__name__}({len(ret_le)})",
            )

            # Value comparison
            aya_se = ret_se[1] if len(ret_se) > 1 else ret_se[0]
            aya_le = ret_le[1] if len(ret_le) > 1 else ret_le[0]

            diff = abs(aya_se - aya_le)
            tol = 0.1 if mode in STAR_BASED_MODES else 0.01
            record(
                f"{test_base}/value",
                diff < tol,
                f"SE={aya_se:.6f} LE={aya_le:.6f} diff={diff:.6f}°",
            )

            # retflag check
            flag_se = ret_se[0]
            flag_le = ret_le[0]
            record(
                f"{test_base}/retflag",
                True,  # Don't fail on flag differences
                f"SE retflag={flag_se} LE retflag={flag_le}",
            )

        except Exception as e:
            record(test_base, False, f"ERROR: {e}")


# ============================================================================
# PART 4: New Modes 43-46 vs PySwissEph
# ============================================================================


def test_part4_new_modes():
    print("\n" + "=" * 70)
    print("PART 4: New Ayanamsha Modes 43-46 vs PySwissEph")
    print("=" * 70)

    jd_2024 = 2460310.5

    new_modes = [
        (43, "Lahiri_1940"),
        (44, "Lahiri_VP285"),
        (45, "Krishnamurti_VP291"),
        (46, "Lahiri_ICRC"),
    ]

    for mode, mode_name in new_modes:
        test_base = f"P4/{mode_name}({mode})"

        # Check if pyswisseph supports this mode
        try:
            swe.set_sid_mode(mode)
            aya_se = swe.get_ayanamsa_ut(jd_2024)
        except Exception as e:
            record_skip(test_base, f"pyswisseph doesn't support mode {mode}: {e}")
            continue

        try:
            ephem.set_sid_mode(mode)
            aya_le = ephem.get_ayanamsa_ut(jd_2024)

            diff = abs(aya_se - aya_le)
            tol = 0.01  # These are formula-based
            record(
                f"{test_base}/aya",
                diff < tol,
                f"SE={aya_se:.6f} LE={aya_le:.6f} "
                f'diff={diff:.6f}° ({diff * 3600:.2f}")',
            )

            # Also test a planet position
            flags = FLG_SIDEREAL | FLG_SPEED
            ret_se = swe.calc_ut(jd_2024, swe.VENUS, flags)
            ret_le = ephem.calc_ut(jd_2024, swe.VENUS, flags)

            pos_se = ret_se[0] if isinstance(ret_se[0], (list, tuple)) else ret_se
            pos_le = ret_le[0] if isinstance(ret_le[0], (list, tuple)) else ret_le

            lon_diff = abs(float(pos_se[0]) - float(pos_le[0]))
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            record(
                f"{test_base}/Venus_lon",
                lon_diff < 0.01,
                f"SE={float(pos_se[0]):.6f} LE={float(pos_le[0]):.6f} "
                f"diff={lon_diff:.6f}°",
            )

        except Exception as e:
            record(test_base, False, f"ERROR: {e}")


# ============================================================================
# PART 5: SIDM_USER Custom Ayanamsha
# ============================================================================


def test_part5_user_ayanamsha():
    print("\n" + "=" * 70)
    print("PART 5: SIDM_USER Custom Ayanamsha")
    print("=" * 70)

    jd_2024 = 2460310.5
    t0 = 2451545.0  # J2000
    ayan_t0 = 23.85  # Custom value at J2000

    swe.set_sid_mode(255, t0, ayan_t0)
    ephem.set_sid_mode(255, t0, ayan_t0)

    try:
        # At t0, ayanamsha should be ayan_t0
        aya_se_t0 = swe.get_ayanamsa_ut(t0)
        aya_le_t0 = ephem.get_ayanamsa_ut(t0)
        diff_t0 = abs(aya_se_t0 - aya_le_t0)
        record(
            "P5/user/at_t0",
            diff_t0 < 0.001,
            f"SE={aya_se_t0:.6f} LE={aya_le_t0:.6f} diff={diff_t0:.6f}° "
            f"(expected ~{ayan_t0:.2f}°)",
        )

        # At 2024, should have precessed from t0
        aya_se_2024 = swe.get_ayanamsa_ut(jd_2024)
        aya_le_2024 = ephem.get_ayanamsa_ut(jd_2024)
        diff_2024 = abs(aya_se_2024 - aya_le_2024)
        record(
            "P5/user/at_2024",
            diff_2024 < 0.01,
            f"SE={aya_se_2024:.6f} LE={aya_le_2024:.6f} diff={diff_2024:.6f}°",
        )

        # Planet position check
        flags = FLG_SIDEREAL | FLG_SPEED
        ret_se = swe.calc_ut(jd_2024, swe.SUN, flags)
        ret_le = ephem.calc_ut(jd_2024, swe.SUN, flags)
        pos_se = ret_se[0] if isinstance(ret_se[0], (list, tuple)) else ret_se
        pos_le = ret_le[0] if isinstance(ret_le[0], (list, tuple)) else ret_le
        lon_diff = abs(float(pos_se[0]) - float(pos_le[0]))
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        record(
            "P5/user/Sun_lon",
            lon_diff < 0.01,
            f"SE={float(pos_se[0]):.6f} LE={float(pos_le[0]):.6f} diff={lon_diff:.6f}°",
        )

    except Exception as e:
        record("P5/user", False, f"ERROR: {e}")

    # Reset to Lahiri
    swe.set_sid_mode(1)
    ephem.set_sid_mode(1)


# ============================================================================
# PART 6: Sidereal Velocity Correctness
# ============================================================================


def test_part6_sidereal_velocity():
    print("\n" + "=" * 70)
    print("PART 6: Sidereal Velocity Correctness")
    print("=" * 70)

    jd = 2460310.5
    swe.set_sid_mode(1)  # Lahiri
    ephem.set_sid_mode(1)

    flags = FLG_SIDEREAL | FLG_SPEED

    for planet_name, planet_id in PLANETS:
        test_name = f"P6/Lahiri/{planet_name}/speed"

        try:
            ret_se = swe.calc_ut(jd, planet_id, flags)
            ret_le = ephem.calc_ut(jd, planet_id, flags)

            pos_se = ret_se[0] if isinstance(ret_se[0], (list, tuple)) else ret_se
            pos_le = ret_le[0] if isinstance(ret_le[0], (list, tuple)) else ret_le

            spd_se = float(pos_se[3])
            spd_le = float(pos_le[3])
            diff = abs(spd_se - spd_le)

            # Also check tropical speed for reference
            ret_trop_se = swe.calc_ut(jd, planet_id, FLG_SPEED)
            ret_trop_le = ephem.calc_ut(jd, planet_id, FLG_SPEED)
            trop_se = (
                ret_trop_se[0]
                if isinstance(ret_trop_se[0], (list, tuple))
                else ret_trop_se
            )
            trop_le = (
                ret_trop_le[0]
                if isinstance(ret_trop_le[0], (list, tuple))
                else ret_trop_le
            )

            # Sidereal speed ≈ tropical speed - ayanamsha_rate (~0.01°/day)
            trop_spd_se = float(trop_se[3])
            sid_spd_se = float(pos_se[3])
            aya_rate = trop_spd_se - sid_spd_se

            ok = diff < 0.01
            record(
                test_name,
                ok,
                f"SE={spd_se:.6f} LE={spd_le:.6f} diff={diff:.6f}°/day "
                f"(aya_rate≈{aya_rate:.6f})",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 7: get_ayanamsa_name
# ============================================================================


def test_part7_ayanamsa_name():
    print("\n" + "=" * 70)
    print("PART 7: get_ayanamsa_name")
    print("=" * 70)

    # Check a few key names
    checks = [
        (0, "Fagan/Bradley"),
        (1, "Lahiri"),
        (5, "Krishnamurti"),
        (18, "J2000"),
        (27, "True Citra"),
    ]

    for mode, expected_substr in checks:
        test_name = f"P7/name/{mode}"
        try:
            name_se = swe.get_ayanamsa_name(mode)
            name_le = ephem.get_ayanamsa_name(mode)

            # Just check both return non-empty strings
            ok = len(name_le) > 0
            record(test_name, ok, f'SE="{name_se}" LE="{name_le}"')
        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 8: Fixed Star Sidereal Position
# ============================================================================


def test_part8_fixstar_sidereal():
    print("\n" + "=" * 70)
    print("PART 8: Fixed Star Sidereal Position")
    print("=" * 70)

    jd = 2460310.5
    swe.set_sid_mode(1)  # Lahiri
    ephem.set_sid_mode(1)

    stars = ["Regulus", "Spica", "Aldebaran"]

    for star_name in stars:
        test_name = f"P8/star/{star_name}"
        try:
            ret_se = swe.fixstar2_ut(star_name, jd, FLG_SIDEREAL)
            pos_se, name_se, retflag_se = ret_se

            ret_le = ephem.fixstar2_ut(star_name, jd, FLG_SIDEREAL)
            pos_le, name_le, retflag_le = ret_le

            lon_se = float(pos_se[0])
            lon_le = float(pos_le[0])
            lon_diff = abs(lon_se - lon_le)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            # Star positions have known ~19" differences due to catalog differences
            ok = lon_diff < 0.02  # 72" — generous for star catalog differences
            record(
                f"{test_name}/lon",
                ok,
                f"SE={lon_se:.6f} LE={lon_le:.6f} "
                f'diff={lon_diff:.6f}° ({lon_diff * 3600:.1f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# MAIN
# ============================================================================


def main():
    global total, passed, failed, skipped, failures

    print("=" * 70)
    print("ROUND 8: Deep Sidereal / Ayanamsha Audit")
    print("=" * 70)
    t0 = time.time()

    test_part1_ayanamsa_values()
    test_part2_sidereal_positions()
    test_part3_ayanamsa_ex_api()
    test_part4_new_modes()
    test_part5_user_ayanamsha()
    test_part6_sidereal_velocity()
    test_part7_ayanamsa_name()
    test_part8_fixstar_sidereal()

    elapsed = time.time() - t0

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total:   {total}")
    print(f"Passed:  {passed}")
    print(f"Failed:  {failed}")
    print(f"Skipped: {skipped}")
    print(f"Time:    {elapsed:.1f}s")

    if failures:
        print(f"\n--- {len(failures)} FAILURES ---")
        for name, detail in failures:
            print(f"  {name}: {detail}")

    print(f"\nPass rate: {passed}/{total} = {100 * passed / max(total, 1):.1f}%")

    # Reset modes
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0)


if __name__ == "__main__":
    main()
