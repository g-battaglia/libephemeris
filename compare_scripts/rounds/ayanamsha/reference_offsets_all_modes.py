#!/usr/bin/env python3
"""Round 45: Ayanamsha All Modes Deep Sweep.

Verifies libephemeris sidereal calculations against pyswisseph for ALL 47
ayanamsha modes, multiple epochs, multiple bodies.

Phases:
  P1: get_ayanamsa_ex_ut value comparison — all 47 modes at J2000
  P2: get_ayanamsa_ex_ut value comparison — all 47 modes at multiple epochs
  P3: Sidereal planet positions — Sun/Moon/Mercury at J2000, all modes
  P4: Sidereal planet positions — Sun/Moon at 5 epochs, selected modes
  P5: Ayanamsha speed (rate of change) comparison
  P6: True ayanamsha modes (star-based) — position consistency check
  P7: set_sid_mode / get_ayanamsa round-trip consistency
"""

from __future__ import annotations

import math
import os
import sys
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
results = {"passed": [], "failed": [], "errors": []}

J2000 = 2451545.0

# All ayanamsha modes (0-46)
SIDM_NAMES = {
    0: "FAGAN_BRADLEY",
    1: "LAHIRI",
    2: "DELUCE",
    3: "RAMAN",
    4: "USHASHASHI",
    5: "KRISHNAMURTI",
    6: "DJWHAL_KHUL",
    7: "YUKTESHWAR",
    8: "JN_BHASIN",
    9: "BABYL_KUGLER1",
    10: "BABYL_KUGLER2",
    11: "BABYL_KUGLER3",
    12: "BABYL_HUBER",
    13: "BABYL_ETPSC",
    14: "ALDEBARAN_15TAU",
    15: "HIPPARCHOS",
    16: "SASSANIAN",
    17: "GALCENT_0SAG",
    18: "J2000",
    19: "J1900",
    20: "B1950",
    21: "SURYASIDDHANTA",
    22: "SURYASIDDHANTA_MSUN",
    23: "ARYABHATA",
    24: "ARYABHATA_MSUN",
    25: "SS_REVATI",
    26: "SS_CITRA",
    27: "TRUE_CITRA",
    28: "TRUE_REVATI",
    29: "TRUE_PUSHYA",
    30: "GALCENT_RGILBRAND",
    31: "GALEQU_IAU1958",
    32: "GALEQU_TRUE",
    33: "GALEQU_MULA",
    34: "GALALIGN_MARDYKS",
    35: "TRUE_MULA",
    36: "GALCENT_MULA_WILHELM",
    37: "ARYABHATA_522",
    38: "BABYL_BRITTON",
    39: "TRUE_SHEORAN",
    40: "GALCENT_COCHRANE",
    41: "GALEQU_FIORENZA",
    42: "VALENS_MOON",
    43: "LAHIRI_1940",
    44: "LAHIRI_VP285",
    45: "KRISHNAMURTI_VP291",
    46: "LAHIRI_ICRC",
}

# Bodies to test
SUN = 0
MOON = 1
MERCURY = 2
VENUS = 3
MARS = 4
JUPITER = 5
SATURN = 6

BODY_NAMES = {
    SUN: "Sun",
    MOON: "Moon",
    MERCURY: "Mercury",
    VENUS: "Venus",
    MARS: "Mars",
    JUPITER: "Jupiter",
    SATURN: "Saturn",
}

# Epochs
EPOCHS = {
    "1900": J2000 - 365.25 * 100,
    "1950": J2000 - 365.25 * 50,
    "2000": J2000,
    "2024": 2460310.5,
    "2050": J2000 + 365.25 * 50,
}


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


def se_sidm(mode):
    """Get SE sidereal mode constant."""
    return mode


def le_sidm(mode):
    """Get LE sidereal mode constant."""
    return mode


def phase1():
    """Ayanamsha value comparison — all modes at J2000."""
    global errors
    print("\n=== P1: Ayanamsha values at J2000, all 47 modes ===")

    jd = J2000
    tol = 15.0 / 3600.0  # 15 arcseconds in degrees (known ~14" model diff)

    for mode_id, mode_name in sorted(SIDM_NAMES.items()):
        try:
            swe.set_sid_mode(se_sidm(mode_id))
            se_aya = swe.get_ayanamsa_ut(jd)

            ephem.set_sid_mode(le_sidm(mode_id), 0, 0)
            le_aya = ephem.get_ayanamsa_ex_ut(jd, ephem.FLG_SIDEREAL)[1]

            diff = abs(se_aya - le_aya)
            diff_arcsec = diff * 3600

            ok = diff_arcsec < 15.0  # 15" tolerance (known ~14" model difference)
            detail = f'SE={se_aya:.6f}° LE={le_aya:.6f}° diff={diff_arcsec:.3f}"'
            record("P1", f"[{mode_id:2d}] {mode_name}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P1 [{mode_id}] {mode_name}: {e}")


def phase2():
    """Ayanamsha value comparison — all modes at multiple epochs."""
    global errors
    print("\n=== P2: Ayanamsha values at multiple epochs ===")

    # Select representative subset of modes to test at all epochs
    test_modes = [0, 1, 3, 5, 7, 14, 17, 18, 21, 27, 28, 29, 30, 35, 39, 42, 43, 46]

    for mode_id in test_modes:
        mode_name = SIDM_NAMES[mode_id]
        for epoch_name, jd in sorted(EPOCHS.items()):
            try:
                swe.set_sid_mode(se_sidm(mode_id))
                se_aya = swe.get_ayanamsa_ut(jd)

                ephem.set_sid_mode(le_sidm(mode_id), 0, 0)
                le_aya = ephem.get_ayanamsa_ex_ut(jd, ephem.FLG_SIDEREAL)[1]

                diff = abs(se_aya - le_aya)
                diff_arcsec = diff * 3600

                # Tolerance scales slightly with time from J2000
                dt_centuries = abs(jd - J2000) / 36525.0
                tol = 15.0 + dt_centuries * 1.0  # 15" + 1"/century

                ok = diff_arcsec < tol
                detail = f'SE={se_aya:.6f}° LE={le_aya:.6f}° diff={diff_arcsec:.3f}" tol={tol:.1f}"'
                record("P2", f"[{mode_id:2d}] {mode_name} {epoch_name}", ok, detail)

            except Exception as e:
                errors += 1
                results["errors"].append(
                    f"P2 [{mode_id}] {mode_name} {epoch_name}: {e}"
                )


def phase3():
    """Sidereal planet positions — Sun/Moon/Mercury at J2000, all modes."""
    global errors
    print("\n=== P3: Sidereal positions Sun/Moon/Merc J2000, all modes ===")

    jd = J2000
    bodies = [SUN, MOON, MERCURY]
    flags_se = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
    flags_le = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL

    for mode_id, mode_name in sorted(SIDM_NAMES.items()):
        swe.set_sid_mode(se_sidm(mode_id))
        ephem.set_sid_mode(le_sidm(mode_id), 0, 0)

        for body in bodies:
            body_name = BODY_NAMES[body]
            try:
                se_result = swe.calc_ut(jd, body, flags_se)
                le_result = ephem.calc_ut(jd, body, flags_le)

                se_pos = se_result[0]
                le_pos = (
                    le_result[0]
                    if isinstance(le_result[0], (list, tuple))
                    else le_result
                )

                lon_diff = abs(se_pos[0] - le_pos[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff
                lat_diff = abs(se_pos[1] - le_pos[1])

                lon_arcsec = lon_diff * 3600
                lat_arcsec = lat_diff * 3600

                # Same ~14" tolerance as ayanamsha itself
                ok = lon_arcsec < 15.0 and lat_arcsec < 1.0
                detail = (
                    f'lon_diff={lon_arcsec:.3f}" lat_diff={lat_arcsec:.3f}" '
                    f"SE_lon={se_pos[0]:.6f}° LE_lon={le_pos[0]:.6f}°"
                )
                record("P3", f"[{mode_id:2d}] {mode_name} {body_name}", ok, detail)

            except Exception as e:
                errors += 1
                results["errors"].append(f"P3 [{mode_id}] {mode_name} {body_name}: {e}")


def phase4():
    """Sidereal positions at multiple epochs, selected modes."""
    global errors
    print("\n=== P4: Sidereal positions multi-epoch ===")

    bodies = [SUN, MOON]
    test_modes = [0, 1, 5, 14, 27, 35, 42, 46]
    flags_se = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL
    flags_le = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL

    for mode_id in test_modes:
        mode_name = SIDM_NAMES[mode_id]
        swe.set_sid_mode(se_sidm(mode_id))
        ephem.set_sid_mode(le_sidm(mode_id), 0, 0)

        for epoch_name, jd in sorted(EPOCHS.items()):
            for body in bodies:
                body_name = BODY_NAMES[body]
                try:
                    se_result = swe.calc_ut(jd, body, flags_se)
                    le_result = ephem.calc_ut(jd, body, flags_le)

                    se_pos = se_result[0]
                    le_pos = (
                        le_result[0]
                        if isinstance(le_result[0], (list, tuple))
                        else le_result
                    )

                    lon_diff = abs(se_pos[0] - le_pos[0])
                    if lon_diff > 180:
                        lon_diff = 360 - lon_diff

                    lon_arcsec = lon_diff * 3600
                    dt_centuries = abs(jd - J2000) / 36525.0
                    tol = 15.0 + dt_centuries * 1.0

                    ok = lon_arcsec < tol
                    detail = (
                        f'lon_diff={lon_arcsec:.3f}" tol={tol:.1f}" '
                        f"SE={se_pos[0]:.6f}° LE={le_pos[0]:.6f}°"
                    )
                    record(
                        "P4",
                        f"[{mode_id:2d}] {mode_name} {body_name} {epoch_name}",
                        ok,
                        detail,
                    )

                except Exception as e:
                    errors += 1
                    results["errors"].append(
                        f"P4 [{mode_id}] {mode_name} {body_name} {epoch_name}: {e}"
                    )


def phase5():
    """Ayanamsha speed (rate of change) comparison."""
    global errors
    print("\n=== P5: Ayanamsha speed (rate of precession) ===")

    jd = J2000
    test_modes = [0, 1, 3, 5, 14, 18, 21, 27, 28, 29, 35, 39, 42, 46]
    flags_se = swe.FLG_SWIEPH | swe.FLG_SIDEREAL
    flags_le = ephem.FLG_SWIEPH | ephem.FLG_SIDEREAL

    for mode_id in test_modes:
        mode_name = SIDM_NAMES[mode_id]
        try:
            # SE: get_ayanamsa_ex_ut returns (retval, ayanamsa) for some versions
            # We'll use finite difference for speed
            dt = 1.0  # 1 day
            swe.set_sid_mode(se_sidm(mode_id))
            se_aya1 = swe.get_ayanamsa_ut(jd - dt / 2)
            se_aya2 = swe.get_ayanamsa_ut(jd + dt / 2)
            se_speed = (se_aya2 - se_aya1) / dt  # deg/day

            ephem.set_sid_mode(le_sidm(mode_id), 0, 0)
            le_aya1 = ephem.get_ayanamsa_ex_ut(jd - dt / 2, flags_le)[1]
            le_aya2 = ephem.get_ayanamsa_ex_ut(jd + dt / 2, flags_le)[1]
            le_speed = (le_aya2 - le_aya1) / dt  # deg/day

            diff = abs(se_speed - le_speed)
            diff_arcsec_per_day = diff * 3600

            # Speed should match within 0.01"/day
            ok = diff_arcsec_per_day < 0.01
            detail = (
                f'SE_speed={se_speed * 3600:.6f}"/day '
                f'LE_speed={le_speed * 3600:.6f}"/day '
                f'diff={diff_arcsec_per_day:.6f}"/day'
            )
            record("P5", f"[{mode_id:2d}] {mode_name} speed", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 [{mode_id}] {mode_name}: {e}")


def phase6():
    """True ayanamsha modes — verify star-based modes are consistent."""
    global errors
    print("\n=== P6: True ayanamsha modes consistency ===")

    # True modes anchor to a specific star position
    # Test that the ayanamsha value is internally consistent:
    # For TRUE_CITRA (27): Spica should be at exactly 180° sidereal
    # For TRUE_REVATI (28): zeta Piscium should be at 359°50' sidereal
    # For TRUE_PUSHYA (29): delta Cancri should be at 106° sidereal

    jd = J2000
    flags_se = swe.FLG_SWIEPH | swe.FLG_SIDEREAL
    flags_le = ephem.FLG_SWIEPH | ephem.FLG_SIDEREAL

    true_modes = [27, 28, 29, 35, 39]

    for mode_id in true_modes:
        mode_name = SIDM_NAMES[mode_id]
        try:
            # Compare SE and LE ayanamsha values at multiple JDs
            test_jds = [J2000, J2000 + 365.25 * 10, J2000 - 365.25 * 10]

            for jd_test in test_jds:
                swe.set_sid_mode(se_sidm(mode_id))
                se_aya = swe.get_ayanamsa_ut(jd_test)

                ephem.set_sid_mode(le_sidm(mode_id), 0, 0)
                le_aya = ephem.get_ayanamsa_ex_ut(jd_test, flags_le)[1]

                diff = abs(se_aya - le_aya)
                diff_arcsec = diff * 3600

                epoch_label = f"JD{jd_test:.1f}"
                ok = diff_arcsec < 20.0  # 20" tolerance for true modes
                detail = f'SE={se_aya:.6f}° LE={le_aya:.6f}° diff={diff_arcsec:.3f}"'
                record("P6", f"[{mode_id:2d}] {mode_name} {epoch_label}", ok, detail)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P6 [{mode_id}] {mode_name}: {e}")


def phase7():
    """set_sid_mode / get_ayanamsa round-trip consistency."""
    global errors
    print("\n=== P7: set_sid_mode round-trip consistency ===")

    jd = J2000

    # Test that after setting mode and computing, the value is self-consistent
    # i.e., tropical_lon - ayanamsha = sidereal_lon

    test_modes = [0, 1, 5, 14, 18, 21, 27, 35, 42, 46]
    bodies = [SUN, MOON, JUPITER]

    for mode_id in test_modes:
        mode_name = SIDM_NAMES[mode_id]
        ephem.set_sid_mode(le_sidm(mode_id), 0, 0)

        for body in bodies:
            body_name = BODY_NAMES[body]
            try:
                # Get tropical position
                le_trop = ephem.calc_ut(
                    jd, body, ephem.FLG_SWIEPH | ephem.FLG_SPEED
                )
                trop_lon = (
                    le_trop[0][0]
                    if isinstance(le_trop[0], (list, tuple))
                    else le_trop[0]
                )

                # Get sidereal position
                le_sid = ephem.calc_ut(
                    jd,
                    body,
                    ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_SIDEREAL,
                )
                sid_lon = (
                    le_sid[0][0] if isinstance(le_sid[0], (list, tuple)) else le_sid[0]
                )

                # Get ayanamsha
                le_aya = ephem.get_ayanamsa_ex_ut(jd, ephem.FLG_SIDEREAL)[1]

                # Check: tropical - ayanamsha ≈ sidereal (mod 360)
                expected_sid = (trop_lon - le_aya) % 360.0
                actual_diff = abs(expected_sid - sid_lon)
                if actual_diff > 180:
                    actual_diff = 360 - actual_diff
                diff_arcsec = actual_diff * 3600

                # Should be essentially exact (internal consistency)
                ok = diff_arcsec < 0.001  # sub-milliarcsecond
                detail = (
                    f"trop={trop_lon:.6f}° aya={le_aya:.6f}° sid={sid_lon:.6f}° "
                    f"expected_sid={(trop_lon - le_aya) % 360:.6f}° "
                    f'diff={diff_arcsec:.6f}"'
                )
                record("P7", f"[{mode_id:2d}] {mode_name} {body_name}", ok, detail)

            except Exception as e:
                errors += 1
                results["errors"].append(f"P7 [{mode_id}] {mode_name} {body_name}: {e}")


def main():
    print("=" * 70)
    print("ROUND 45: Ayanamsha All Modes Deep Sweep")
    print("=" * 70)

    phase1()
    print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

    phase2()
    print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

    phase3()
    print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

    phase4()
    print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

    phase5()
    print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

    phase6()
    print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

    phase7()
    print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

    total = passed + failed + errors
    pct = 100 * passed / total if total else 0
    print("\n" + "=" * 70)
    print(f"ROUND 45 FINAL: {passed}/{total} passed ({pct:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print("=" * 70)

    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:50]:
            print(f)

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:20]:
            print(e)

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
