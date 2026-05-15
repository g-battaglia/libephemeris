#!/usr/bin/env python3
"""Round 105: Sidereal Time at Geographic Locations

Tests sidtime() and sidtime0() across epochs and locations.
Compares GMST and LST values between SE and LE.
"""

from __future__ import annotations
import sys, os, time
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

EPOCHS = [2451545.0 + i * 365.25 for i in range(50)]  # 50 years
LONGITUDES = [
    0.0,
    30.0,
    60.0,
    90.0,
    120.0,
    150.0,
    180.0,
    -30.0,
    -60.0,
    -90.0,
    -120.0,
    -150.0,
]

TOL_SECONDS = 0.1  # 0.1 second of time


def run_tests():
    passed = 0
    failed = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 105: Sidereal Time at Geographic Locations")
    print("=" * 80)

    # PART 1: sidtime (GMST)
    print("\n--- PART 1: sidtime (GMST) ---")
    p1_pass = 0
    p1_fail = 0

    for jd in EPOCHS:
        total += 1
        se_st = swe.sidtime(jd)
        le_st = ephem.sidtime(jd)
        diff_sec = abs(se_st - le_st) * 3600  # hours -> seconds

        if diff_sec <= TOL_SECONDS:
            p1_pass += 1
            passed += 1
        else:
            p1_fail += 1
            failed += 1
            if len(fail_details) < 5:
                fail_details.append(
                    f"  FAIL [sidtime] JD={jd:.1f}: SE={se_st:.10f}h LE={le_st:.10f}h diff={diff_sec:.4f}s"
                )

    print(f"  Part 1: {p1_pass}/{p1_pass + p1_fail} passed")

    # PART 2: sidtime0 (with obliquity and nutation)
    print("\n--- PART 2: sidtime0 (with eps/dpsi) ---")
    p2_pass = 0
    p2_fail = 0

    for jd in EPOCHS[:20]:
        total += 1
        # Get obliquity and nutation from SE
        nut = swe.calc_ut(jd, -1, 2)  # FLG_SWIEPH
        eps = nut[0][0]  # true obliquity
        dpsi = nut[0][2]  # nutation in longitude (degrees)

        se_st0 = swe.sidtime0(jd, eps, dpsi)
        le_st0 = ephem.sidtime0(jd, eps, dpsi)
        diff_sec = abs(se_st0 - le_st0) * 3600

        if diff_sec <= TOL_SECONDS:
            p2_pass += 1
            passed += 1
        else:
            p2_fail += 1
            failed += 1
            if len(fail_details) < 10:
                fail_details.append(
                    f"  FAIL [sidtime0] JD={jd:.1f}: SE={se_st0:.10f}h LE={le_st0:.10f}h diff={diff_sec:.4f}s"
                )

    print(f"  Part 2: {p2_pass}/{p2_pass + p2_fail} passed")

    # PART 3: Local sidereal time (GMST + longitude offset)
    print("\n--- PART 3: Local Sidereal Time ---")
    p3_pass = 0
    p3_fail = 0

    for jd in EPOCHS[:20]:
        se_gmst = swe.sidtime(jd)
        le_gmst = ephem.sidtime(jd)

        for lon in LONGITUDES:
            total += 1
            se_lst = (se_gmst + lon / 15.0) % 24.0
            le_lst = (le_gmst + lon / 15.0) % 24.0
            diff_sec = abs(se_lst - le_lst) * 3600

            if diff_sec > 12 * 3600:  # wrap around
                diff_sec = 24 * 3600 - diff_sec

            if diff_sec <= TOL_SECONDS:
                p3_pass += 1
                passed += 1
            else:
                p3_fail += 1
                failed += 1

    print(f"  Part 3: {p3_pass}/{p3_pass + p3_fail} passed")

    # PART 4: Consistency with house calculation ARMC
    print("\n--- PART 4: Sidereal Time vs Houses ARMC ---")
    p4_pass = 0
    p4_fail = 0

    for jd in EPOCHS[:20]:
        for lon in [0.0, 45.0, 90.0, -45.0, -90.0]:
            total += 1
            lat = 45.0

            try:
                se_cusps, se_ascmc = swe.houses(jd, lat, lon, b"P")
                le_result = ephem.houses(jd, lat, lon, ord("P"))
                le_ascmc = le_result[1]

                # ARMC is ascmc[2]
                diff = abs(se_ascmc[2] - le_ascmc[2])
                if diff > 180:
                    diff = 360 - diff

                if diff < 0.01:  # 0.01 degree = 36"
                    p4_pass += 1
                    passed += 1
                else:
                    p4_fail += 1
                    failed += 1
                    if len(fail_details) < 15:
                        fail_details.append(
                            f"  FAIL [ARMC] JD={jd:.1f} lon={lon}: SE={se_ascmc[2]:.6f} LE={le_ascmc[2]:.6f}"
                        )
            except Exception:
                p4_fail += 1
                failed += 1

    print(f"  Part 4: {p4_pass}/{p4_pass + p4_fail} passed")

    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details:
            print(d)

    pct = 100 * passed / max(1, total)
    print()
    print("=" * 80)
    print(f"ROUND 105 FINAL: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
