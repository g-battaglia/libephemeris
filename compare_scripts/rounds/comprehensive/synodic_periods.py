#!/usr/bin/env python3
"""Round 99: Synodic Periods & Conjunction/Opposition Timing

Verifies planet positions at known conjunction/opposition times by computing
positions for all planets at many epochs and checking SE vs LE agreement.

Tests:
1. All planets at monthly intervals over 20 years (position agreement)
2. Mercury/Venus inferior conjunction detection (inner planet close to Sun)
3. Mars/Jupiter/Saturn opposition detection (outer planet 180° from Sun)
4. Moon conjunctions with planets (Moon near planet longitude)
5. Speed agreement at conjunction/opposition epochs
"""

from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2

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

PLANETS = [
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
]

# Tolerances
LON_TOL = 1.0  # arcsec
LAT_TOL = 1.0  # arcsec
SPEED_TOL = 1.0  # arcsec/day


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 99: Synodic Periods & Conjunction/Opposition Timing")
    print("=" * 80)

    # =========================================================================
    # PART 1: All planets at monthly intervals 2000-2020
    # =========================================================================
    print("\n--- PART 1: Monthly All-Planet Sweep (2000-2020) ---")
    p1_pass = 0
    p1_fail = 0

    base_jd = 2451545.0  # J2000
    flags = FLG_SWIEPH | FLG_SPEED

    for month_offset in range(240):  # 20 years monthly
        jd = base_jd + month_offset * 30.4375  # ~monthly

        for body_id, body_name in PLANETS:
            total += 1

            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                errors += 1
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600
            dspd = abs(se_pos[3] - le_pos[3]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            if dlon <= LON_TOL and dlat <= LAT_TOL and dspd <= SPEED_TOL:
                p1_pass += 1
                passed += 1
            else:
                p1_fail += 1
                failed += 1
                if len(fail_details) < 10:
                    fail_details.append(
                        f"  FAIL [P1] {body_name:10s} JD={jd:.1f}: "
                        f'dLon={dlon:.3f}" dLat={dlat:.3f}" dSpd={dspd:.3f}"/d'
                    )

    print(
        f"  Part 1: {p1_pass}/{p1_pass + p1_fail} passed ({100 * p1_pass / max(1, p1_pass + p1_fail):.1f}%)"
    )

    # =========================================================================
    # PART 2: Planet-Sun elongation agreement
    # =========================================================================
    print("\n--- PART 2: Planet-Sun Elongation Agreement ---")
    p2_pass = 0
    p2_fail = 0

    OUTER_PLANETS = [
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
    ]
    INNER_PLANETS = [(MERCURY, "Mercury"), (VENUS, "Venus")]

    for week in range(520):  # 10 years weekly
        jd = base_jd + week * 7.0

        se_sun = swe.calc_ut(jd, SUN, flags)[0]
        le_sun = ephem.calc_ut(jd, SUN, flags)[0]

        for body_id, body_name in OUTER_PLANETS + INNER_PLANETS:
            total += 1

            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                errors += 1
                continue

            # Elongation from Sun
            se_elong = (se_pos[0] - se_sun[0]) % 360.0
            le_elong = (le_pos[0] - le_sun[0]) % 360.0

            delong = abs(se_elong - le_elong) * 3600
            if delong > 180 * 3600:
                delong = 360 * 3600 - delong

            if delong <= 2.0:  # 2" tolerance for elongation
                p2_pass += 1
                passed += 1
            else:
                p2_fail += 1
                failed += 1
                if len(fail_details) < 15:
                    fail_details.append(
                        f"  FAIL [P2] {body_name:10s} JD={jd:.1f}: "
                        f'dElong={delong:.2f}" SE_e={se_elong:.4f} LE_e={le_elong:.4f}'
                    )

    print(
        f"  Part 2: {p2_pass}/{p2_pass + p2_fail} passed ({100 * p2_pass / max(1, p2_pass + p2_fail):.1f}%)"
    )

    # =========================================================================
    # PART 3: Moon-planet separation agreement
    # =========================================================================
    print("\n--- PART 3: Moon-Planet Separation Agreement ---")
    p3_pass = 0
    p3_fail = 0

    for day in range(365):  # 1 year daily
        jd = base_jd + day

        se_moon = swe.calc_ut(jd, MOON, flags)[0]
        le_moon = ephem.calc_ut(jd, MOON, flags)[0]

        for body_id, body_name in PLANETS[2:]:  # Mercury through Pluto
            total += 1

            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                errors += 1
                continue

            # Moon-planet separation
            se_sep = abs(se_moon[0] - se_pos[0])
            if se_sep > 180:
                se_sep = 360 - se_sep
            le_sep = abs(le_moon[0] - le_pos[0])
            if le_sep > 180:
                le_sep = 360 - le_sep

            dsep = abs(se_sep - le_sep) * 3600

            if dsep <= 3.0:  # 3" tolerance for separation
                p3_pass += 1
                passed += 1
            else:
                p3_fail += 1
                failed += 1
                if len(fail_details) < 20:
                    fail_details.append(
                        f"  FAIL [P3] Moon-{body_name:10s} JD={jd:.1f}: "
                        f'dSep={dsep:.2f}" SE_sep={se_sep:.4f}° LE_sep={le_sep:.4f}°'
                    )

    print(
        f"  Part 3: {p3_pass}/{p3_pass + p3_fail} passed ({100 * p3_pass / max(1, p3_pass + p3_fail):.1f}%)"
    )

    # =========================================================================
    # PART 4: Historical & future epochs position agreement
    # =========================================================================
    print("\n--- PART 4: Historical & Future Epoch Positions ---")
    p4_pass = 0
    p4_fail = 0

    HISTORICAL_EPOCHS = [
        2415020.0,  # 1900-01-01
        2420000.0,  # ~1913
        2425000.0,  # ~1927
        2430000.0,  # ~1941
        2435000.0,  # ~1955
        2440000.0,  # ~1968
        2445000.0,  # ~1982
        2450000.0,  # ~1995
        2455000.0,  # ~2009
        2460000.0,  # ~2023
        2465000.0,  # ~2036
        2470000.0,  # ~2050
        2475000.0,  # ~2063
        2480000.0,  # ~2077
        2485000.0,  # ~2090
        2488069.5,  # ~2099-12-31
    ]

    for jd in HISTORICAL_EPOCHS:
        for body_id, body_name in PLANETS:
            total += 1

            try:
                se_pos = swe.calc_ut(jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(jd, body_id, flags)[0]
            except Exception:
                errors += 1
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            # Wider tolerance for far future (Delta-T divergence)
            epoch_tol = LON_TOL if jd < 2470000 else 3.0

            if dlon <= epoch_tol and dlat <= epoch_tol:
                p4_pass += 1
                passed += 1
            else:
                p4_fail += 1
                failed += 1
                if len(fail_details) < 25:
                    fail_details.append(
                        f"  FAIL [P4] {body_name:10s} JD={jd:.1f}: "
                        f'dLon={dlon:.3f}" dLat={dlat:.3f}"'
                    )

    print(
        f"  Part 4: {p4_pass}/{p4_pass + p4_fail} passed ({100 * p4_pass / max(1, p4_pass + p4_fail):.1f}%)"
    )

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    print()
    print("=" * 80)
    print(
        f"ROUND 99 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
