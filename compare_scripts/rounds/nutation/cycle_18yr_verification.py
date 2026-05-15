#!/usr/bin/env python3
"""Round 106: Nutation 18.6-Year Cycle Verification

Verifies nutation angles (dpsi, deps) and their effect on planet positions
across a full 18.6-year lunar nodal regression cycle.
"""

from __future__ import annotations
import sys, os, time, math
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2
FLG_NONUT = 64

# Sample every 30 days over 20 years (covers one full 18.6-year cycle)
EPOCHS = [2451545.0 + i * 30.0 for i in range(243)]

# Tolerances
NUT_TOL = 0.5  # arcsec for nutation angles
POS_TOL = 1.0  # arcsec for planet positions
NUT_EFFECT_TOL = 2.0  # arcsec for nutation effect on positions


def run_tests():
    passed = 0
    failed = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 106: Nutation 18.6-Year Cycle Verification")
    print("=" * 80)

    # =========================================================================
    # PART 1: Nutation angles (dpsi, deps) from ECL_NUT
    # =========================================================================
    print("\n--- PART 1: Nutation Angles (dpsi, deps) ---")
    p1_pass = 0
    p1_fail = 0
    worst_dpsi = 0
    worst_deps = 0

    for jd in EPOCHS:
        total += 1
        se_nut = swe.calc_ut(jd, -1, FLG_SWIEPH)[0]  # ECL_NUT
        le_nut = ephem.calc_ut(jd, -1, FLG_SWIEPH)[0]

        # se_nut: (true_obl, mean_obl, dpsi, deps, 0, 0)
        ddpsi = abs(se_nut[2] - le_nut[2]) * 3600  # arcsec
        ddeps = abs(se_nut[3] - le_nut[3]) * 3600

        worst_dpsi = max(worst_dpsi, ddpsi)
        worst_deps = max(worst_deps, ddeps)

        if ddpsi <= NUT_TOL and ddeps <= NUT_TOL:
            p1_pass += 1
            passed += 1
        else:
            p1_fail += 1
            failed += 1
            if len(fail_details) < 5:
                fail_details.append(
                    f'  FAIL [NUT] JD={jd:.1f}: ddpsi={ddpsi:.3f}" ddeps={ddeps:.3f}"'
                )

    print(
        f"  Part 1: {p1_pass}/{p1_pass + p1_fail} passed "
        f'(worst dpsi={worst_dpsi:.3f}" deps={worst_deps:.3f}")'
    )

    # =========================================================================
    # PART 2: True vs mean obliquity agreement
    # =========================================================================
    print("\n--- PART 2: True & Mean Obliquity ---")
    p2_pass = 0
    p2_fail = 0
    worst_obl = 0

    for jd in EPOCHS:
        total += 1
        se_nut = swe.calc_ut(jd, -1, FLG_SWIEPH)[0]
        le_nut = ephem.calc_ut(jd, -1, FLG_SWIEPH)[0]

        d_true = abs(se_nut[0] - le_nut[0]) * 3600
        d_mean = abs(se_nut[1] - le_nut[1]) * 3600
        worst_obl = max(worst_obl, d_true, d_mean)

        if d_true <= NUT_TOL and d_mean <= NUT_TOL:
            p2_pass += 1
            passed += 1
        else:
            p2_fail += 1
            failed += 1

    print(f'  Part 2: {p2_pass}/{p2_pass + p2_fail} passed (worst={worst_obl:.3f}")')

    # =========================================================================
    # PART 3: Nutation effect on Sun position (with vs without NONUT)
    # =========================================================================
    print("\n--- PART 3: Nutation Effect on Sun Position ---")
    p3_pass = 0
    p3_fail = 0

    for jd in EPOCHS:
        total += 1
        flags_nut = FLG_SWIEPH | FLG_SPEED
        flags_nonut = FLG_SWIEPH | FLG_SPEED | FLG_NONUT

        se_nut_pos = swe.calc_ut(jd, 0, flags_nut)[0]
        se_nonut_pos = swe.calc_ut(jd, 0, flags_nonut)[0]
        le_nut_pos = ephem.calc_ut(jd, 0, flags_nut)[0]
        le_nonut_pos = ephem.calc_ut(jd, 0, flags_nonut)[0]

        # Nutation effect in longitude
        se_nut_effect = (se_nut_pos[0] - se_nonut_pos[0]) * 3600
        le_nut_effect = (le_nut_pos[0] - le_nonut_pos[0]) * 3600

        diff = abs(se_nut_effect - le_nut_effect)

        if diff <= NUT_EFFECT_TOL:
            p3_pass += 1
            passed += 1
        else:
            p3_fail += 1
            failed += 1
            if len(fail_details) < 10:
                fail_details.append(
                    f'  FAIL [SUN_NUT] JD={jd:.1f}: SE_eff={se_nut_effect:.3f}" '
                    f'LE_eff={le_nut_effect:.3f}" diff={diff:.3f}"'
                )

    print(f"  Part 3: {p3_pass}/{p3_pass + p3_fail} passed")

    # =========================================================================
    # PART 4: Moon position with and without NONUT
    # =========================================================================
    print("\n--- PART 4: Nutation Effect on Moon Position ---")
    p4_pass = 0
    p4_fail = 0

    for jd in EPOCHS:
        total += 1
        flags_nut = FLG_SWIEPH | FLG_SPEED
        flags_nonut = FLG_SWIEPH | FLG_SPEED | FLG_NONUT

        se_nut_pos = swe.calc_ut(jd, 1, flags_nut)[0]
        se_nonut_pos = swe.calc_ut(jd, 1, flags_nonut)[0]
        le_nut_pos = ephem.calc_ut(jd, 1, flags_nut)[0]
        le_nonut_pos = ephem.calc_ut(jd, 1, flags_nonut)[0]

        se_nut_effect = (se_nut_pos[0] - se_nonut_pos[0]) * 3600
        le_nut_effect = (le_nut_pos[0] - le_nonut_pos[0]) * 3600

        diff = abs(se_nut_effect - le_nut_effect)

        if diff <= NUT_EFFECT_TOL:
            p4_pass += 1
            passed += 1
        else:
            p4_fail += 1
            failed += 1
            if len(fail_details) < 15:
                fail_details.append(
                    f'  FAIL [MOON_NUT] JD={jd:.1f}: SE={se_nut_effect:.3f}" '
                    f'LE={le_nut_effect:.3f}" diff={diff:.3f}"'
                )

    print(f"  Part 4: {p4_pass}/{p4_pass + p4_fail} passed")

    # =========================================================================
    # PART 5: Nutation amplitude range check
    # =========================================================================
    print("\n--- PART 5: Nutation Amplitude Range Check ---")
    p5_pass = 0
    p5_fail = 0
    le_dpsi_min = 999
    le_dpsi_max = -999

    for jd in EPOCHS:
        total += 1
        le_nut = ephem.calc_ut(jd, -1, FLG_SWIEPH)[0]
        dpsi = le_nut[2] * 3600  # arcsec

        le_dpsi_min = min(le_dpsi_min, dpsi)
        le_dpsi_max = max(le_dpsi_max, dpsi)

        # IAU 2000B nutation: dpsi should be roughly -17" to +17"
        if -20.0 <= dpsi <= 20.0:
            p5_pass += 1
            passed += 1
        else:
            p5_fail += 1
            failed += 1

    print(
        f"  Part 5: {p5_pass}/{p5_pass + p5_fail} passed "
        f'(dpsi range: {le_dpsi_min:.2f}" to {le_dpsi_max:.2f}")'
    )

    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    pct = 100 * passed / max(1, total)
    print()
    print("=" * 80)
    print(f"ROUND 106 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
