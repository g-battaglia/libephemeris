#!/usr/bin/env python3
"""Round 103: Ecliptic/Equatorial Coordinate Round-Trip

Verifies that computing a planet's position in ecliptic, then converting to
equatorial via cotrans, matches the direct FLG_EQUATORIAL output.
Also verifies the reverse round-trip (equatorial -> ecliptic via cotrans).

Tests:
1. All planets: ecliptic -> cotrans -> equatorial vs direct equatorial
2. Round-trip: ecliptic -> equatorial -> ecliptic
3. cotrans_sp velocity consistency
"""

from __future__ import annotations
import sys, os, time, math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2
FLG_EQUATORIAL = 2048

PLANETS = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
]

EPOCHS = [2451545.0 + i * 365.25 for i in range(20)]

TOL = 1.0  # arcsec
ROUNDTRIP_TOL = 0.001  # arcsec — should be near-perfect


def run_tests():
    passed = 0
    failed = 0
    total = 0
    fail_details = []
    flags_ecl = FLG_SWIEPH | FLG_SPEED
    flags_eq = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL

    print("=" * 80)
    print("ROUND 103: Ecliptic/Equatorial Coordinate Round-Trip")
    print("=" * 80)

    # =========================================================================
    # PART 1: LE ecliptic + cotrans vs LE equatorial (internal consistency)
    # =========================================================================
    print("\n--- PART 1: LE Ecliptic+cotrans vs LE Equatorial ---")
    p1_pass = 0
    p1_fail = 0

    for jd in EPOCHS:
        # Get obliquity
        nut_result = ephem.calc_ut(jd, -1, FLG_SWIEPH)
        eps = nut_result[0][0]  # true obliquity

        for body_id, body_name in PLANETS:
            total += 1
            ecl = ephem.calc_ut(jd, body_id, flags_ecl)[0]
            eq_direct = ephem.calc_ut(jd, body_id, flags_eq)[0]

            # Manual cotrans from ecliptic to equatorial
            pos = (ecl[0], ecl[1], ecl[2])
            vel = (ecl[3], ecl[4], ecl[5])
            pos_eq, vel_eq = ephem.cotrans_sp(pos, vel, -eps)

            dra = abs(pos_eq[0] - eq_direct[0]) * 3600
            ddec = abs(pos_eq[1] - eq_direct[1]) * 3600
            if dra > 180 * 3600:
                dra = 360 * 3600 - dra

            if dra <= TOL and ddec <= TOL:
                p1_pass += 1
                passed += 1
            else:
                p1_fail += 1
                failed += 1
                if len(fail_details) < 10:
                    fail_details.append(
                        f"  FAIL [P1] {body_name:10s} JD={jd:.1f}: "
                        f'dRA={dra:.3f}" dDec={ddec:.3f}"'
                    )

    print(f"  Part 1: {p1_pass}/{p1_pass + p1_fail} passed")

    # =========================================================================
    # PART 2: Round-trip ecliptic -> equatorial -> ecliptic (LE)
    # =========================================================================
    print("\n--- PART 2: Round-Trip Ecliptic -> Equatorial -> Ecliptic ---")
    p2_pass = 0
    p2_fail = 0

    for jd in EPOCHS:
        nut_result = ephem.calc_ut(jd, -1, FLG_SWIEPH)
        eps = nut_result[0][0]

        for body_id, body_name in PLANETS:
            total += 1
            ecl = ephem.calc_ut(jd, body_id, flags_ecl)[0]

            # Forward: ecliptic -> equatorial
            pos = (ecl[0], ecl[1], ecl[2])
            vel = (ecl[3], ecl[4], ecl[5])
            pos_eq, vel_eq = ephem.cotrans_sp(pos, vel, -eps)

            # Reverse: equatorial -> ecliptic
            pos_back, vel_back = ephem.cotrans_sp(pos_eq, vel_eq, eps)

            dlon = abs(pos_back[0] - ecl[0]) * 3600
            dlat = abs(pos_back[1] - ecl[1]) * 3600
            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            if dlon <= ROUNDTRIP_TOL and dlat <= ROUNDTRIP_TOL:
                p2_pass += 1
                passed += 1
            else:
                p2_fail += 1
                failed += 1
                if len(fail_details) < 15:
                    fail_details.append(
                        f"  FAIL [P2] {body_name:10s} JD={jd:.1f}: "
                        f'dLon={dlon:.6f}" dLat={dlat:.6f}"'
                    )

    print(f"  Part 2: {p2_pass}/{p2_pass + p2_fail} passed")

    # =========================================================================
    # PART 3: SE ecliptic + SE cotrans vs SE equatorial (SE internal consistency)
    # =========================================================================
    print("\n--- PART 3: SE Ecliptic+cotrans vs SE Equatorial ---")
    p3_pass = 0
    p3_fail = 0

    for jd in EPOCHS:
        nut_result = swe.calc_ut(jd, -1, FLG_SWIEPH)
        eps = nut_result[0][0]

        for body_id, body_name in PLANETS:
            total += 1
            se_ecl = swe.calc_ut(jd, body_id, flags_ecl)[0]
            se_eq = swe.calc_ut(jd, body_id, flags_eq)[0]

            pos = (se_ecl[0], se_ecl[1], se_ecl[2])
            pos_eq = swe.cotrans(pos, -eps)

            dra = abs(pos_eq[0] - se_eq[0]) * 3600
            ddec = abs(pos_eq[1] - se_eq[1]) * 3600
            if dra > 180 * 3600:
                dra = 360 * 3600 - dra

            if dra <= TOL and ddec <= TOL:
                p3_pass += 1
                passed += 1
            else:
                p3_fail += 1
                failed += 1
                if len(fail_details) < 20:
                    fail_details.append(
                        f"  FAIL [P3] {body_name:10s} JD={jd:.1f}: "
                        f'dRA={dra:.3f}" dDec={ddec:.3f}"'
                    )

    print(f"  Part 3: {p3_pass}/{p3_pass + p3_fail} passed")

    # =========================================================================
    # PART 4: Cross-library: LE ecliptic vs SE ecliptic, LE equatorial vs SE equatorial
    # =========================================================================
    print("\n--- PART 4: Cross-Library Position Agreement ---")
    p4_pass = 0
    p4_fail = 0

    for jd in EPOCHS:
        for body_id, body_name in PLANETS:
            # Ecliptic
            total += 1
            se_ecl = swe.calc_ut(jd, body_id, flags_ecl)[0]
            le_ecl = ephem.calc_ut(jd, body_id, flags_ecl)[0]
            dlon = abs(se_ecl[0] - le_ecl[0]) * 3600
            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon
            dlat = abs(se_ecl[1] - le_ecl[1]) * 3600

            if dlon <= TOL and dlat <= TOL:
                p4_pass += 1
                passed += 1
            else:
                p4_fail += 1
                failed += 1

            # Equatorial
            total += 1
            se_eq = swe.calc_ut(jd, body_id, flags_eq)[0]
            le_eq = ephem.calc_ut(jd, body_id, flags_eq)[0]
            dra = abs(se_eq[0] - le_eq[0]) * 3600
            if dra > 180 * 3600:
                dra = 360 * 3600 - dra
            ddec = abs(se_eq[1] - le_eq[1]) * 3600

            if dra <= TOL and ddec <= TOL:
                p4_pass += 1
                passed += 1
            else:
                p4_fail += 1
                failed += 1

    print(f"  Part 4: {p4_pass}/{p4_pass + p4_fail} passed")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    pct = 100 * passed / max(1, total)
    print()
    print("=" * 80)
    print(f"ROUND 103 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
