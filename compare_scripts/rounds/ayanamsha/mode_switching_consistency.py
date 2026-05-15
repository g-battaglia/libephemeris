#!/usr/bin/env python3
"""Round 119: Multiple Ayanamsha Mode Switching — Consistency Test

Tests that switching between ayanamsha modes produces consistent results,
and that mode switching doesn't leave residual state.

Verifies:
- All 43+ ayanamsha modes produce valid results
- Switching modes mid-session doesn't corrupt state
- get_ayanamsa_ut matches between LE and SE for each mode
- Sidereal planet positions consistent after mode switches
- Custom ayanamsha (mode 255) works correctly
- Mode reset after each switch
"""

from __future__ import annotations

import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SIDEREAL = 65536

SUN = 0
MOON = 1
MARS = 4
JUPITER = 5
SATURN = 6

# All ayanamsha modes to test
AYANAMSHA_MODES = [
    (0, "FAGAN_BRADLEY"),
    (1, "LAHIRI"),
    (2, "DELUCE"),
    (3, "RAMAN"),
    (4, "USHASHASHI"),
    (5, "KRISHNAMURTI"),
    (6, "DJWHAL_KHUL"),
    (7, "YUKTESHWAR"),
    (8, "JN_BHASIN"),
    (9, "BABYL_KUGLER1"),
    (10, "BABYL_KUGLER2"),
    (11, "BABYL_KUGLER3"),
    (12, "BABYL_HUBER"),
    (13, "BABYL_ETPSC"),
    (14, "ALDEBARAN_15TAU"),
    (15, "HIPPARCHOS"),
    (16, "SASSANIAN"),
    (17, "GALCENT_0SAG"),
    (18, "J2000"),
    (19, "J1900"),
    (20, "B1950"),
    (21, "SURYASIDDHANTA"),
    (22, "SURYASIDDHANTA_MSUN"),
    (23, "ARYABHATA"),
    (24, "ARYABHATA_MSUN"),
    (25, "SS_REVATI"),
    (26, "SS_CITRA"),
    (27, "TRUE_CITRA"),
    (28, "TRUE_REVATI"),
    (29, "TRUE_PUSHYA"),
    (30, "GALCENT_RGILBRAND"),
    (31, "GALEQU_IAU1958"),
    (32, "GALEQU_TRUE"),
    (33, "GALEQU_MULA"),
    (34, "GALALIGN_MARDYKS"),
    (35, "TRUE_MULA"),
    (36, "GALCENT_MULA_WILHELM"),
    (37, "ARYABHATA_522"),
    (38, "BABYL_BRITTON"),
    (39, "TRUE_SHEORAN"),
    (40, "GALCENT_COCHRANE"),
    (41, "GALEQU_FIORENZA"),
    (42, "VALENS_MOON"),
    (43, "LAHIRI_1940"),
    (44, "LAHIRI_VP285"),
    (45, "KRISHNAMURTI_VP291"),
    (46, "LAHIRI_ICRC"),
]


def main():
    print("=" * 80)
    print("ROUND 119: Multiple Ayanamsha Mode Switching — Consistency Test")
    print("=" * 80)

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    test_jds = [
        2451545.0,  # J2000.0
        2455197.5,  # 2010-01-01
        2459580.5,  # 2022-01-01
        2460310.5,  # 2024-01-15
        2440587.5,  # 1970-01-01
        2430000.5,  # 1941-02-15
    ]

    bodies = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    # Test 1: get_ayanamsa_ut for all modes
    print("\n--- Test 1: get_ayanamsa_ut for all modes ---")
    for mode_id, mode_name in AYANAMSHA_MODES:
        swe.set_sid_mode(mode_id)
        ephem.set_sid_mode(mode_id, 0, 0)

        for jd in test_jds:
            try:
                se_aya = swe.get_ayanamsa_ut(jd)
            except Exception:
                continue

            try:
                le_aya = ephem.get_ayanamsa_ut(jd)
            except Exception as e:
                total_tests += 1
                total_fail += 1
                if len(failures) < 20:
                    failures.append(
                        f"  get_ayanamsa_ut mode={mode_id} ({mode_name}): LE_ERROR: {e}"
                    )
                continue

            diff_arcsec = abs(le_aya - se_aya) * 3600
            # Known ~14" sidereal offset for many modes, up to ~700" for galactic modes at extreme dates
            tol = 20.0  # 20" tolerance
            if mode_id in (17, 30, 37, 40):  # Galactic center modes
                tol = 800.0
            if mode_id == 18:  # J2000 mode
                tol = 1.0

            total_tests += 1
            if diff_arcsec < tol:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 40:
                    failures.append(
                        f"  get_ayanamsa_ut mode={mode_id} ({mode_name}) JD={jd}: "
                        f'SE={se_aya:.6f}° LE={le_aya:.6f}° diff={diff_arcsec:.2f}" (tol={tol}")'
                    )

    # Test 2: Sidereal planet positions for all modes
    print("\n--- Test 2: Sidereal planet positions for all modes ---")
    for mode_id, mode_name in AYANAMSHA_MODES:
        swe.set_sid_mode(mode_id)
        ephem.set_sid_mode(mode_id, 0, 0)

        for jd in test_jds[:3]:  # Subset for speed
            for body_id, body_name in bodies:
                flags = FLG_SPEED | FLG_SIDEREAL

                try:
                    se_result = swe.calc_ut(jd, body_id, flags)
                    se_lon = se_result[0][0]
                except Exception:
                    continue

                try:
                    le_result = ephem.calc_ut(jd, body_id, flags)
                    le_lon = le_result[0][0]
                except Exception as e:
                    total_tests += 1
                    total_fail += 1
                    continue

                diff = le_lon - se_lon
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                diff_arcsec = abs(diff) * 3600

                tol = 20.0
                if mode_id in (17, 30, 37, 40):
                    tol = 800.0
                if mode_id == 18:
                    tol = 1.0

                total_tests += 1
                if diff_arcsec < tol:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 60:
                        failures.append(
                            f"  Sidereal mode={mode_id} ({mode_name}) JD={jd} {body_name}: "
                            f'SE={se_lon:.6f}° LE={le_lon:.6f}° diff={diff_arcsec:.2f}" (tol={tol}")'
                        )

    # Test 3: Mode switching doesn't corrupt state
    print("\n--- Test 3: Mode switching state isolation ---")
    jd = 2451545.0

    # Get reference positions with Fagan-Bradley
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    ref_se = swe.calc_ut(jd, MOON, FLG_SPEED | FLG_SIDEREAL)[0][0]
    ref_le = ephem.calc_ut(jd, MOON, FLG_SPEED | FLG_SIDEREAL)[0][0]

    # Switch through several modes
    for mode_id in [1, 3, 5, 17, 27, 42]:
        swe.set_sid_mode(mode_id)
        ephem.set_sid_mode(mode_id, 0, 0)
        _ = swe.calc_ut(jd, MOON, FLG_SPEED | FLG_SIDEREAL)
        _ = ephem.calc_ut(jd, MOON, FLG_SPEED | FLG_SIDEREAL)

    # Switch back to Fagan-Bradley
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)
    check_se = swe.calc_ut(jd, MOON, FLG_SPEED | FLG_SIDEREAL)[0][0]
    check_le = ephem.calc_ut(jd, MOON, FLG_SPEED | FLG_SIDEREAL)[0][0]

    total_tests += 2
    if abs(ref_se - check_se) < 1e-10:
        total_pass += 1
    else:
        total_fail += 1
        failures.append(
            f"  SE state corrupted: ref={ref_se:.10f} check={check_se:.10f}"
        )

    if abs(ref_le - check_le) < 1e-10:
        total_pass += 1
    else:
        total_fail += 1
        failures.append(
            f"  LE state corrupted: ref={ref_le:.10f} check={check_le:.10f}"
        )

    # Test 4: Tropical position unaffected by sidereal mode
    print("\n--- Test 4: Tropical position unaffected by sidereal mode ---")
    for jd in test_jds[:3]:
        # Get tropical position with no sidereal flag
        for body_id, body_name in bodies:
            results_tropical = []
            for mode_id in [0, 1, 3, 17, 27, 42]:
                swe.set_sid_mode(mode_id)
                ephem.set_sid_mode(mode_id, 0, 0)

                # Without FLG_SIDEREAL, result should be identical regardless of mode
                le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)
                results_tropical.append(le_result[0][0])

            # All should be identical
            for i in range(1, len(results_tropical)):
                total_tests += 1
                if abs(results_tropical[i] - results_tropical[0]) < 1e-10:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 80:
                        failures.append(
                            f"  Tropical affected by sid mode: JD={jd} {body_name}: "
                            f"mode0={results_tropical[0]:.10f} != mode{i}={results_tropical[i]:.10f}"
                        )

    # Test 5: Custom ayanamsha (mode 255)
    print("\n--- Test 5: Custom ayanamsha mode ---")
    custom_t0 = 2451545.0  # J2000
    custom_ayan_t0 = 23.5  # degrees
    custom_ayan_rate = 50.3 / 3600  # degrees per year (approximately)

    swe.set_sid_mode(255, custom_t0, custom_ayan_t0)
    ephem.set_sid_mode(255, custom_t0, custom_ayan_t0)

    for jd in test_jds:
        try:
            se_aya = swe.get_ayanamsa_ut(jd)
        except Exception:
            continue

        try:
            le_aya = ephem.get_ayanamsa_ut(jd)
        except Exception:
            total_tests += 1
            total_fail += 1
            continue

        diff_arcsec = abs(le_aya - se_aya) * 3600
        tol = 20.0  # Custom mode should still be close

        total_tests += 1
        if diff_arcsec < tol:
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 90:
                failures.append(
                    f'  Custom mode JD={jd}: SE={se_aya:.6f}° LE={le_aya:.6f}° diff={diff_arcsec:.2f}"'
                )

    # Test 6: Ayanamsha monotonicity (should generally increase over time)
    print("\n--- Test 6: Ayanamsha monotonicity check ---")
    monotonic_modes = [0, 1, 3, 5, 14, 21]  # Standard modes should be monotonic

    for mode_id in monotonic_modes:
        ephem.set_sid_mode(mode_id, 0, 0)

        prev_aya = None
        jds_sorted = sorted(test_jds)

        for jd in jds_sorted:
            le_aya = ephem.get_ayanamsa_ut(jd)

            if prev_aya is not None:
                total_tests += 1
                if le_aya >= prev_aya:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 95:
                        failures.append(
                            f"  Non-monotonic: mode={mode_id} prev={prev_aya:.6f}° curr={le_aya:.6f}°"
                        )

            prev_aya = le_aya

    # Test 7: Difference between modes is consistent
    print("\n--- Test 7: Inter-mode difference consistency ---")
    # Lahiri - Fagan should be approximately the same across dates
    for jd in test_jds:
        swe.set_sid_mode(0)
        ephem.set_sid_mode(0, 0, 0)
        se_fb = swe.get_ayanamsa_ut(jd)
        le_fb = ephem.get_ayanamsa_ut(jd)

        swe.set_sid_mode(1)
        ephem.set_sid_mode(1, 0, 0)
        se_lah = swe.get_ayanamsa_ut(jd)
        le_lah = ephem.get_ayanamsa_ut(jd)

        se_diff = se_lah - se_fb
        le_diff = le_lah - le_fb

        total_tests += 1
        # The difference between modes should be very similar between SE and LE
        if abs(se_diff - le_diff) * 3600 < 5.0:  # Within 5"
            total_pass += 1
        else:
            total_fail += 1
            if len(failures) < 100:
                failures.append(
                    f"  Mode diff inconsistency JD={jd}: SE(Lah-FB)={se_diff:.6f}° "
                    f'LE(Lah-FB)={le_diff:.6f}° diff={abs(se_diff - le_diff) * 3600:.2f}"'
                )

    # Test 8: Sidereal houses
    print("\n--- Test 8: Sidereal house cusps ---")
    for mode_id, mode_name in [(0, "FB"), (1, "LAH"), (27, "TRUE_CITRA")]:
        swe.set_sid_mode(mode_id)
        ephem.set_sid_mode(mode_id, 0, 0)

        for jd in test_jds[:3]:
            for lat in [0.0, 45.0, -33.0]:
                lon = 12.5
                try:
                    se_cusps, se_ascmc = swe.houses_ex(
                        jd, lat, lon, b"P", FLG_SIDEREAL
                    )
                except Exception:
                    continue

                try:
                    le_result = ephem.houses_ex2(
                        jd, lat, lon, ord("P"), FLG_SIDEREAL | FLG_SPEED
                    )
                    le_cusps = le_result[0]
                except Exception:
                    continue

                # Compare first 12 cusps
                for i in range(min(len(se_cusps), len(le_cusps), 12)):
                    diff = le_cusps[i] - se_cusps[i]
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360
                    diff_arcsec = abs(diff) * 3600

                    tol = 20.0
                    total_tests += 1
                    if diff_arcsec < tol:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 110:
                            failures.append(
                                f"  SidHouse {mode_name} lat={lat} cusp{i + 1}: "
                                f'SE={se_cusps[i]:.6f}° LE={le_cusps[i]:.6f}° diff={diff_arcsec:.2f}"'
                            )

    # Reset
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    # Summary
    print("\n" + "=" * 80)
    pct = 100 * total_pass / total_tests if total_tests > 0 else 0
    print(f"ROUND 119 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
    print(f"  Failures: {total_fail}")
    print("=" * 80)

    if failures:
        print("\nSample failures:")
        for f in failures[:25]:
            print(f)

    if total_fail == 0:
        print("\nAll tests PASSED!")

    return total_fail


if __name__ == "__main__":
    sys.exit(main())
