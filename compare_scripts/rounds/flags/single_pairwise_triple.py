#!/usr/bin/env python3
"""
Round 10: Combined Flags Stress Test
======================================
Tests many flag combinations to find edge cases:
  P1: All single-flag modes for Sun/Moon/Jupiter
  P2: Pairwise flag combinations (EQUATORIAL+J2000, SIDEREAL+SPEED, etc.)
  P3: Triple flag combinations
  P4: Flag preservation in retflag
  P5: FLG_TRUEPOS / FLG_NOABERR / FLG_NOGDEFL
  P6: Heliocentric + Barycentric with various sub-flags
  P7: Distance/radius values across flag modes
  P8: FLG_RADIANS output conversion
"""

from __future__ import annotations

import math
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


# ============================================================================
# CONSTANTS
# ============================================================================

JD = 2460310.5  # 2024-Jan-15

# Flag definitions — use pyswisseph constants to avoid errors
FLG_SPEED = swe.FLG_SPEED
FLG_TOPOCTR = swe.FLG_TOPOCTR
FLG_EQUATORIAL = swe.FLG_EQUATORIAL
FLG_SIDEREAL = swe.FLG_SIDEREAL
FLG_J2000 = swe.FLG_J2000
FLG_NONUT = swe.FLG_NONUT
FLG_TRUEPOS = swe.FLG_TRUEPOS
FLG_NOABERR = swe.FLG_NOABERR
FLG_NOGDEFL = swe.FLG_NOGDEFL
FLG_HELCTR = swe.FLG_HELCTR
FLG_BARYCTR = swe.FLG_BARYCTR
FLG_XYZ = swe.FLG_XYZ
FLG_RADIANS = swe.FLG_RADIANS

BODIES = [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.JUPITER, "Jupiter"),
]


# ============================================================================
# HELPERS
# ============================================================================


def compare_pos(
    test_name,
    pos_se,
    pos_le,
    lon_tol=0.0005,
    lat_tol=0.0005,
    speed_tol=0.001,
    is_xyz=False,
    is_radians=False,
):
    """Compare SE and LE position tuples, record results."""
    labels = (
        ("x", "y", "z", "vx", "vy", "vz")
        if is_xyz
        else ("lon", "lat", "dist", "speed_lon", "speed_lat", "speed_dist")
    )

    for i in range(min(6, len(pos_se), len(pos_le))):
        val_se = float(pos_se[i])
        val_le = float(pos_le[i])
        diff = abs(val_se - val_le)

        # For longitude, handle wrap-around
        if i == 0 and not is_xyz and not is_radians:
            if diff > 180:
                diff = 360 - diff

        if i == 0 and is_radians and not is_xyz:
            if diff > math.pi:
                diff = 2 * math.pi - diff

        # Determine tolerance
        if i < 2:
            tol = lon_tol
        elif i == 2:
            tol = 0.001  # distance AU
        else:
            tol = speed_tol

        ok = diff < tol
        record(
            f"{test_name}/{labels[i]}",
            ok,
            f"SE={val_se:.8f} LE={val_le:.8f} diff={diff:.8f}",
        )


# ============================================================================
# PART 1: All single-flag modes
# ============================================================================


def test_part1_single_flags():
    print("\n" + "=" * 70)
    print("PART 1: Single Flag Modes")
    print("=" * 70)

    single_flags = [
        ("bare", 0),
        ("SPEED", FLG_SPEED),
        ("EQUATORIAL", FLG_EQUATORIAL),
        ("J2000", FLG_J2000),
        ("NONUT", FLG_NONUT),
        ("TRUEPOS", FLG_TRUEPOS),
        ("NOABERR", FLG_NOABERR),
        ("NOGDEFL", FLG_NOGDEFL),
        ("HELCTR", FLG_HELCTR),
        ("BARYCTR", FLG_BARYCTR),
        ("XYZ", FLG_XYZ),
    ]

    for body_id, body_name in BODIES:
        for flag_name, flag_val in single_flags:
            # Skip helio Sun (meaningless)
            if body_name == "Sun" and flag_name == "HELCTR":
                continue
            test_name = f"P1/{body_name}/{flag_name}"
            try:
                pos_se, rf_se = swe.calc_ut(JD, body_id, flag_val)
                pos_le, rf_le = ephem.calc_ut(JD, body_id, flag_val)

                is_xyz = bool(flag_val & FLG_XYZ)
                lon_se = float(pos_se[0])
                lon_le = float(pos_le[0])

                if is_xyz:
                    diff = abs(lon_se - lon_le)
                else:
                    diff = abs(lon_se - lon_le)
                    if diff > 180:
                        diff = 360 - diff

                tol = 0.001 if body_name == "Moon" else 0.0005
                ok = diff < tol
                record(
                    test_name,
                    ok,
                    f'SE[0]={lon_se:.6f} LE[0]={lon_le:.6f} diff={diff:.6f}° ({diff * 3600:.2f}")',
                )

            except Exception as e:
                record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 2: Pairwise flag combinations
# ============================================================================


def test_part2_pairwise():
    print("\n" + "=" * 70)
    print("PART 2: Pairwise Flag Combinations")
    print("=" * 70)

    pairs = [
        ("SPEED+EQ", FLG_SPEED | FLG_EQUATORIAL),
        ("SPEED+J2000", FLG_SPEED | FLG_J2000),
        ("SPEED+NONUT", FLG_SPEED | FLG_NONUT),
        ("SPEED+TRUEPOS", FLG_SPEED | FLG_TRUEPOS),
        ("SPEED+NOABERR", FLG_SPEED | FLG_NOABERR),
        ("SPEED+NOGDEFL", FLG_SPEED | FLG_NOGDEFL),
        ("EQ+J2000", FLG_EQUATORIAL | FLG_J2000),
        ("EQ+NONUT", FLG_EQUATORIAL | FLG_NONUT),
        ("J2000+NOABERR", FLG_J2000 | FLG_NOABERR),
        ("NONUT+NOABERR", FLG_NONUT | FLG_NOABERR),
        ("TRUEPOS+EQ", FLG_TRUEPOS | FLG_EQUATORIAL),
        ("NOABERR+NOGDEFL", FLG_NOABERR | FLG_NOGDEFL),
        ("HELCTR+SPEED", FLG_HELCTR | FLG_SPEED),
        ("BARYCTR+SPEED", FLG_BARYCTR | FLG_SPEED),
        ("HELCTR+EQ", FLG_HELCTR | FLG_EQUATORIAL),
        ("BARYCTR+EQ", FLG_BARYCTR | FLG_EQUATORIAL),
        ("HELCTR+J2000", FLG_HELCTR | FLG_J2000),
        ("BARYCTR+J2000", FLG_BARYCTR | FLG_J2000),
    ]

    body_id, body_name = swe.JUPITER, "Jupiter"

    for pair_name, flag_val in pairs:
        test_name = f"P2/{body_name}/{pair_name}"
        try:
            pos_se, _ = swe.calc_ut(JD, body_id, flag_val)
            pos_le, _ = ephem.calc_ut(JD, body_id, flag_val)

            is_xyz = bool(flag_val & FLG_XYZ)
            val_se = float(pos_se[0])
            val_le = float(pos_le[0])

            if is_xyz:
                diff = abs(val_se - val_le)
            else:
                diff = abs(val_se - val_le)
                if diff > 180:
                    diff = 360 - diff

            ok = diff < 0.0005
            record(
                test_name,
                ok,
                f'SE[0]={val_se:.6f} LE[0]={val_le:.6f} diff={diff:.6f}° ({diff * 3600:.2f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # Also test Moon for a subset
    body_id, body_name = swe.MOON, "Moon"
    moon_pairs = [
        ("SPEED+EQ", FLG_SPEED | FLG_EQUATORIAL),
        ("SPEED+J2000", FLG_SPEED | FLG_J2000),
        ("EQ+J2000", FLG_EQUATORIAL | FLG_J2000),
        ("TRUEPOS+EQ", FLG_TRUEPOS | FLG_EQUATORIAL),
        ("NOABERR+EQ", FLG_NOABERR | FLG_EQUATORIAL),
    ]

    for pair_name, flag_val in moon_pairs:
        test_name = f"P2/{body_name}/{pair_name}"
        try:
            pos_se, _ = swe.calc_ut(JD, body_id, flag_val)
            pos_le, _ = ephem.calc_ut(JD, body_id, flag_val)

            val_se = float(pos_se[0])
            val_le = float(pos_le[0])
            diff = abs(val_se - val_le)
            if diff > 180:
                diff = 360 - diff

            ok = diff < 0.001
            record(
                test_name,
                ok,
                f'SE[0]={val_se:.6f} LE[0]={val_le:.6f} diff={diff:.6f}° ({diff * 3600:.2f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 3: Triple flag combinations
# ============================================================================


def test_part3_triple():
    print("\n" + "=" * 70)
    print("PART 3: Triple Flag Combinations")
    print("=" * 70)

    triples = [
        ("SPEED+EQ+J2000", FLG_SPEED | FLG_EQUATORIAL | FLG_J2000),
        ("SPEED+EQ+NONUT", FLG_SPEED | FLG_EQUATORIAL | FLG_NONUT),
        ("SPEED+J2000+NOABERR", FLG_SPEED | FLG_J2000 | FLG_NOABERR),
        ("SPEED+TRUEPOS+EQ", FLG_SPEED | FLG_TRUEPOS | FLG_EQUATORIAL),
        ("SPEED+NOABERR+NOGDEFL", FLG_SPEED | FLG_NOABERR | FLG_NOGDEFL),
        ("SPEED+EQ+NOABERR", FLG_SPEED | FLG_EQUATORIAL | FLG_NOABERR),
        ("HELCTR+SPEED+EQ", FLG_HELCTR | FLG_SPEED | FLG_EQUATORIAL),
        ("HELCTR+SPEED+J2000", FLG_HELCTR | FLG_SPEED | FLG_J2000),
        ("BARYCTR+SPEED+EQ", FLG_BARYCTR | FLG_SPEED | FLG_EQUATORIAL),
        ("SPEED+J2000+NOGDEFL", FLG_SPEED | FLG_J2000 | FLG_NOGDEFL),
    ]

    for body_id, body_name in [(swe.MOON, "Moon"), (swe.MARS, "Mars")]:
        for triple_name, flag_val in triples:
            test_name = f"P3/{body_name}/{triple_name}"
            try:
                pos_se, _ = swe.calc_ut(JD, body_id, flag_val)
                pos_le, _ = ephem.calc_ut(JD, body_id, flag_val)

                val_se = float(pos_se[0])
                val_le = float(pos_le[0])
                diff = abs(val_se - val_le)
                if diff > 180:
                    diff = 360 - diff

                tol = 0.001 if body_name == "Moon" else 0.0005
                ok = diff < tol
                record(
                    test_name,
                    ok,
                    f'SE[0]={val_se:.6f} LE[0]={val_le:.6f} diff={diff:.6f}° ({diff * 3600:.2f}")',
                )

            except Exception as e:
                record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 4: Flag preservation in retflag
# ============================================================================


def test_part4_retflag():
    print("\n" + "=" * 70)
    print("PART 4: Flag Preservation in retflag")
    print("=" * 70)

    # Flags that should be preserved in retflag
    flag_tests = [
        ("SPEED", FLG_SPEED),
        ("EQUATORIAL", FLG_EQUATORIAL),
        ("J2000", FLG_J2000),
        ("NONUT", FLG_NONUT),
        ("TRUEPOS", FLG_TRUEPOS),
        ("NOABERR", FLG_NOABERR),
        ("NOGDEFL", FLG_NOGDEFL),
        ("XYZ", FLG_XYZ),
        ("RADIANS", FLG_RADIANS),
        ("HELCTR", FLG_HELCTR),
        ("BARYCTR", FLG_BARYCTR),
    ]

    body_id = swe.MARS

    for flag_name, flag_val in flag_tests:
        test_name = f"P4/retflag/{flag_name}"
        try:
            _, rf_le = ephem.calc_ut(JD, body_id, flag_val)

            # Check that the flag bit is set in retflag
            flag_preserved = bool(rf_le & flag_val)
            record(
                test_name,
                flag_preserved,
                f"iflag=0x{flag_val:06x} retflag=0x{rf_le:06x} "
                f"bit_set={flag_preserved}",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # Combined flags
    combined = FLG_SPEED | FLG_EQUATORIAL | FLG_J2000
    test_name = "P4/retflag/SPEED+EQ+J2000"
    try:
        _, rf_le = ephem.calc_ut(JD, body_id, combined)
        all_set = (rf_le & combined) == combined
        record(
            test_name,
            all_set,
            f"iflag=0x{combined:06x} retflag=0x{rf_le:06x} all_bits={all_set}",
        )
    except Exception as e:
        record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 5: TRUEPOS / NOABERR / NOGDEFL effect on positions
# ============================================================================


def test_part5_aberration_deflection():
    print("\n" + "=" * 70)
    print("PART 5: Aberration & Deflection Flag Effects")
    print("=" * 70)

    body_id, body_name = swe.MARS, "Mars"

    # Get baseline (apparent position — includes aberration + deflection)
    pos_apparent_se, _ = swe.calc_ut(JD, body_id, FLG_SPEED)
    pos_apparent_le, _ = ephem.calc_ut(JD, body_id, FLG_SPEED)

    # Get TRUEPOS (geometric — no aberration, no deflection)
    pos_true_se, _ = swe.calc_ut(JD, body_id, FLG_SPEED | FLG_TRUEPOS)
    pos_true_le, _ = ephem.calc_ut(JD, body_id, FLG_SPEED | FLG_TRUEPOS)

    # Get NOABERR (no aberration, but WITH deflection)
    pos_noaberr_se, _ = swe.calc_ut(JD, body_id, FLG_SPEED | FLG_NOABERR)
    pos_noaberr_le, _ = ephem.calc_ut(JD, body_id, FLG_SPEED | FLG_NOABERR)

    # Get NOGDEFL (WITH aberration, no deflection)
    pos_nogdefl_se, _ = swe.calc_ut(JD, body_id, FLG_SPEED | FLG_NOGDEFL)
    pos_nogdefl_le, _ = ephem.calc_ut(JD, body_id, FLG_SPEED | FLG_NOGDEFL)

    # 1. Verify apparent positions match
    diff = abs(float(pos_apparent_se[0]) - float(pos_apparent_le[0]))
    if diff > 180:
        diff = 360 - diff
    record("P5/apparent/Mars", diff < 0.0005, f'diff={diff:.6f}° ({diff * 3600:.2f}")')

    # 2. Verify TRUEPOS positions match
    diff = abs(float(pos_true_se[0]) - float(pos_true_le[0]))
    if diff > 180:
        diff = 360 - diff
    record("P5/truepos/Mars", diff < 0.0005, f'diff={diff:.6f}° ({diff * 3600:.2f}")')

    # 3. Verify NOABERR positions match
    diff = abs(float(pos_noaberr_se[0]) - float(pos_noaberr_le[0]))
    if diff > 180:
        diff = 360 - diff
    record("P5/noaberr/Mars", diff < 0.0005, f'diff={diff:.6f}° ({diff * 3600:.2f}")')

    # 4. Verify NOGDEFL positions match
    diff = abs(float(pos_nogdefl_se[0]) - float(pos_nogdefl_le[0]))
    if diff > 180:
        diff = 360 - diff
    record("P5/nogdefl/Mars", diff < 0.0005, f'diff={diff:.6f}° ({diff * 3600:.2f}")')

    # 5. Aberration effect should be ~20" (annual aberration)
    aberr_se = abs(float(pos_apparent_se[0]) - float(pos_noaberr_se[0]))
    if aberr_se > 180:
        aberr_se = 360 - aberr_se
    aberr_le = abs(float(pos_apparent_le[0]) - float(pos_noaberr_le[0]))
    if aberr_le > 180:
        aberr_le = 360 - aberr_le
    aberr_diff = abs(aberr_se - aberr_le)
    record(
        "P5/aberr_effect/Mars",
        aberr_diff < 0.001,
        f'SE_aberr={aberr_se * 3600:.2f}" LE_aberr={aberr_le * 3600:.2f}" diff={aberr_diff * 3600:.2f}"',
    )

    # 6. Deflection effect (gravitational light bending near Sun, usually <1")
    defl_se = abs(float(pos_apparent_se[0]) - float(pos_nogdefl_se[0]))
    if defl_se > 180:
        defl_se = 360 - defl_se
    defl_le = abs(float(pos_apparent_le[0]) - float(pos_nogdefl_le[0]))
    if defl_le > 180:
        defl_le = 360 - defl_le
    defl_diff = abs(defl_se - defl_le)
    record(
        "P5/defl_effect/Mars",
        defl_diff < 0.001,
        f'SE_defl={defl_se * 3600:.2f}" LE_defl={defl_le * 3600:.2f}" diff={defl_diff * 3600:.2f}"',
    )

    # Repeat for Sun (aberration should be ~0 for Sun, but TRUEPOS differs)
    for body_id2, bname2 in [(swe.SUN, "Sun"), (swe.MOON, "Moon")]:
        for mode, mflag in [
            ("truepos", FLG_TRUEPOS),
            ("noaberr", FLG_NOABERR),
            ("nogdefl", FLG_NOGDEFL),
        ]:
            try:
                pos_se2, _ = swe.calc_ut(JD, body_id2, FLG_SPEED | mflag)
                pos_le2, _ = ephem.calc_ut(JD, body_id2, FLG_SPEED | mflag)
                diff2 = abs(float(pos_se2[0]) - float(pos_le2[0]))
                if diff2 > 180:
                    diff2 = 360 - diff2
                tol = 0.001 if bname2 == "Moon" else 0.0005
                record(
                    f"P5/{mode}/{bname2}",
                    diff2 < tol,
                    f'diff={diff2:.6f}° ({diff2 * 3600:.2f}")',
                )
            except Exception as e:
                record(f"P5/{mode}/{bname2}", False, f"ERROR: {e}")


# ============================================================================
# PART 6: Heliocentric + Barycentric with sub-flags
# ============================================================================


def test_part6_helio_bary():
    print("\n" + "=" * 70)
    print("PART 6: Heliocentric & Barycentric with Sub-flags")
    print("=" * 70)

    bodies_hb = [
        (swe.MOON, "Moon"),
        (swe.MARS, "Mars"),
        (swe.JUPITER, "Jupiter"),
        (swe.SATURN, "Saturn"),
    ]

    for center_name, center_flag in [("helio", FLG_HELCTR), ("bary", FLG_BARYCTR)]:
        for body_id, body_name in bodies_hb:
            for sub_name, sub_flag in [
                ("bare", 0),
                ("SPEED", FLG_SPEED),
                ("EQ", FLG_EQUATORIAL),
                ("J2000", FLG_J2000),
                ("SPEED+EQ", FLG_SPEED | FLG_EQUATORIAL),
            ]:
                flags = center_flag | sub_flag
                test_name = f"P6/{center_name}/{body_name}/{sub_name}"
                try:
                    pos_se, _ = swe.calc_ut(JD, body_id, flags)
                    pos_le, _ = ephem.calc_ut(JD, body_id, flags)

                    val_se = float(pos_se[0])
                    val_le = float(pos_le[0])
                    diff = abs(val_se - val_le)
                    if not (flags & FLG_XYZ) and diff > 180:
                        diff = 360 - diff

                    tol = 0.001
                    ok = diff < tol
                    record(
                        test_name,
                        ok,
                        f'SE={val_se:.6f} LE={val_le:.6f} diff={diff:.6f}° ({diff * 3600:.2f}")',
                    )

                except Exception as e:
                    record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 7: Distance consistency across flag modes
# ============================================================================


def test_part7_distances():
    print("\n" + "=" * 70)
    print("PART 7: Distance Consistency Across Flag Modes")
    print("=" * 70)

    body_id, body_name = swe.MARS, "Mars"

    modes = [
        ("default", FLG_SPEED),
        ("J2000", FLG_SPEED | FLG_J2000),
        ("NONUT", FLG_SPEED | FLG_NONUT),
        ("EQUATORIAL", FLG_SPEED | FLG_EQUATORIAL),
        ("EQ+J2000", FLG_SPEED | FLG_EQUATORIAL | FLG_J2000),
        ("TRUEPOS", FLG_SPEED | FLG_TRUEPOS),
        ("NOABERR", FLG_SPEED | FLG_NOABERR),
    ]

    for mode_name, flags in modes:
        test_name = f"P7/dist/{body_name}/{mode_name}"
        try:
            pos_se, _ = swe.calc_ut(JD, body_id, flags)
            pos_le, _ = ephem.calc_ut(JD, body_id, flags)

            dist_se = float(pos_se[2])
            dist_le = float(pos_le[2])
            diff = abs(dist_se - dist_le)

            ok = diff < 0.0001  # AU
            record(
                test_name,
                ok,
                f"SE={dist_se:.8f} LE={dist_le:.8f} diff={diff:.8f} AU",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 8: FLG_RADIANS output
# ============================================================================


def test_part8_radians():
    print("\n" + "=" * 70)
    print("PART 8: FLG_RADIANS Output Conversion")
    print("=" * 70)

    for body_id, body_name in BODIES:
        test_name = f"P8/radians/{body_name}"
        try:
            # Get degrees
            pos_deg_se, _ = swe.calc_ut(JD, body_id, FLG_SPEED)
            pos_deg_le, _ = ephem.calc_ut(JD, body_id, FLG_SPEED)

            # Get radians
            pos_rad_se, _ = swe.calc_ut(JD, body_id, FLG_SPEED | FLG_RADIANS)
            pos_rad_le, _ = ephem.calc_ut(JD, body_id, FLG_SPEED | FLG_RADIANS)

            # Verify SE radians match
            lon_rad_se = float(pos_rad_se[0])
            lon_rad_le = float(pos_rad_le[0])
            diff_rad = abs(lon_rad_se - lon_rad_le)
            if diff_rad > math.pi:
                diff_rad = 2 * math.pi - diff_rad

            tol = 0.00002 if body_name == "Moon" else 0.00001  # ~1" in radians
            ok = diff_rad < tol
            record(
                f"{test_name}/rad_match",
                ok,
                f"SE={lon_rad_se:.8f} LE={lon_rad_le:.8f} diff={diff_rad:.8f} rad",
            )

            # Verify LE radians are consistent with LE degrees
            lon_deg_le = float(pos_deg_le[0])
            expected_rad = math.radians(lon_deg_le)
            conv_diff = abs(lon_rad_le - expected_rad)
            if conv_diff > math.pi:
                conv_diff = 2 * math.pi - conv_diff

            ok_conv = conv_diff < 1e-10
            record(
                f"{test_name}/deg_rad_consistent",
                ok_conv,
                f"deg={lon_deg_le:.6f} rad={lon_rad_le:.8f} expected_rad={expected_rad:.8f} diff={conv_diff:.2e}",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # XYZ + RADIANS
    for body_id, body_name in [(swe.MARS, "Mars")]:
        test_name = f"P8/xyz_radians/{body_name}"
        try:
            pos_se, _ = swe.calc_ut(
                JD, body_id, FLG_SPEED | FLG_XYZ | FLG_RADIANS
            )
            pos_le, _ = ephem.calc_ut(
                JD, body_id, FLG_SPEED | FLG_XYZ | FLG_RADIANS
            )

            # XYZ positions should not change with RADIANS (already Cartesian)
            for i, label in enumerate(["x", "y", "z"]):
                val_se = float(pos_se[i])
                val_le = float(pos_le[i])
                diff = abs(val_se - val_le)
                ok = diff < 0.001
                record(
                    f"{test_name}/{label}",
                    ok,
                    f"SE={val_se:.8f} LE={val_le:.8f} diff={diff:.8f}",
                )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# MAIN
# ============================================================================


def main():
    global total, passed, failed, skipped, failures

    print("=" * 70)
    print("ROUND 10: Combined Flags Stress Test")
    print("=" * 70)
    t0 = time.time()

    test_part1_single_flags()
    test_part2_pairwise()
    test_part3_triple()
    test_part4_retflag()
    test_part5_aberration_deflection()
    test_part6_helio_bary()
    test_part7_distances()
    test_part8_radians()

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

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
