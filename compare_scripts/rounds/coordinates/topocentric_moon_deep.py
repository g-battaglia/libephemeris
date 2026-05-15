#!/usr/bin/env python3
"""Round 97: Topocentric Moon Parallax — Deep Sweep

The Moon's parallax is the largest topocentric correction (~1°). This round
verifies topocentric Moon positions at many geographic locations and epochs.

Tests:
1. Multiple geographic locations (cities worldwide)
2. Different altitudes (sea level, mountain, aircraft)
3. Multiple epochs across a synodic month
4. Comparison of lon, lat, dist, and speeds
5. Both ecliptic and equatorial topocentric output
"""

from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SWIEPH = 2
FLG_TOPOCTR = 32768
FLG_EQUATORIAL = 2048
MOON = 1

# Geographic locations: (lat, lon, alt_m, name)
LOCATIONS = [
    (51.5074, -0.1278, 11, "London"),
    (40.7128, -74.0060, 10, "New York"),
    (35.6762, 139.6503, 40, "Tokyo"),
    (-33.8688, 151.2093, 58, "Sydney"),
    (-22.9068, -43.1729, 11, "Rio"),
    (55.7558, 37.6173, 156, "Moscow"),
    (28.6139, 77.2090, 216, "Delhi"),
    (1.3521, 103.8198, 15, "Singapore"),
    (64.1466, -21.9426, 0, "Reykjavik"),
    (-54.8019, -68.3030, 18, "Ushuaia"),
    (0.0, 0.0, 0, "Null Island"),
    (71.0, 25.0, 0, "Hammerfest"),
    (-77.85, 166.67, 0, "McMurdo"),
    (27.9881, 86.9250, 8849, "Everest"),  # High altitude
    (30.0, 35.0, -430, "Dead Sea"),  # Below sea level
]

# Epochs: one synodic month around J2000
BASE_JD = 2451545.0
EPOCHS = [BASE_JD + i * 2.5 for i in range(12)]  # ~30 days, every 2.5 days
# Add some other interesting epochs
EPOCHS += [
    2440000.0,  # 1968
    2460000.0,  # 2023
    2470000.0,  # 2050
]

# Tolerances
LON_TOL = 2.0  # arcsec (topocentric models may differ slightly)
LAT_TOL = 2.0  # arcsec
DIST_TOL = 5e-5  # AU (~7.5 km)
SPEED_TOL = 2.0  # arcsec/day
RA_TOL = 2.0  # arcsec
DEC_TOL = 2.0  # arcsec


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 97: Topocentric Moon Parallax — Deep Sweep")
    print("=" * 80)

    # =========================================================================
    # PART 1: Topocentric Moon ecliptic positions
    # =========================================================================
    print("\n--- PART 1: Topocentric Moon Ecliptic Positions ---")
    p1_pass = 0
    p1_fail = 0
    p1_err = 0

    for loc_lat, loc_lon, loc_alt, loc_name in LOCATIONS:
        swe.set_topo(loc_lon, loc_lat, loc_alt)
        ephem.set_topo(loc_lon, loc_lat, loc_alt)

        for jd in EPOCHS:
            total += 1
            flags = FLG_SWIEPH | FLG_SPEED | FLG_TOPOCTR

            try:
                se_result = swe.calc_ut(jd, MOON, flags)
                se_pos = se_result[0]
            except Exception as e:
                p1_err += 1
                errors += 1
                continue

            try:
                le_result = ephem.calc_ut(jd, MOON, flags)
                le_pos = le_result[0]
            except Exception as e:
                p1_err += 1
                errors += 1
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600
            ddist = abs(se_pos[2] - le_pos[2])
            dspd = abs(se_pos[3] - le_pos[3]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            ok = (
                dlon <= LON_TOL
                and dlat <= LAT_TOL
                and ddist <= DIST_TOL
                and dspd <= SPEED_TOL
            )

            if ok:
                p1_pass += 1
                passed += 1
            else:
                p1_fail += 1
                failed += 1
                if len(fail_details) < 15:
                    fail_details.append(
                        f"  FAIL [ECL] {loc_name:12s} JD={jd:.1f}: "
                        f'dLon={dlon:.2f}" dLat={dlat:.2f}" dDist={ddist:.6f} '
                        f'dSpd={dspd:.2f}"/d'
                    )

    print(
        f"  Part 1: {p1_pass}/{p1_pass + p1_fail + p1_err} passed, {p1_fail} fail, {p1_err} err"
    )

    # =========================================================================
    # PART 2: Topocentric Moon equatorial positions
    # =========================================================================
    print("\n--- PART 2: Topocentric Moon Equatorial Positions ---")
    p2_pass = 0
    p2_fail = 0
    p2_err = 0

    for loc_lat, loc_lon, loc_alt, loc_name in LOCATIONS:
        swe.set_topo(loc_lon, loc_lat, loc_alt)
        ephem.set_topo(loc_lon, loc_lat, loc_alt)

        for jd in EPOCHS[:8]:  # Subset for equatorial
            total += 1
            flags = FLG_SWIEPH | FLG_SPEED | FLG_TOPOCTR | FLG_EQUATORIAL

            try:
                se_result = swe.calc_ut(jd, MOON, flags)
                se_pos = se_result[0]
            except Exception as e:
                p2_err += 1
                errors += 1
                continue

            try:
                le_result = ephem.calc_ut(jd, MOON, flags)
                le_pos = le_result[0]
            except Exception as e:
                p2_err += 1
                errors += 1
                continue

            dra = abs(se_pos[0] - le_pos[0]) * 3600
            ddec = abs(se_pos[1] - le_pos[1]) * 3600

            if dra > 180 * 3600:
                dra = 360 * 3600 - dra

            ok = dra <= RA_TOL and ddec <= DEC_TOL

            if ok:
                p2_pass += 1
                passed += 1
            else:
                p2_fail += 1
                failed += 1
                if len(fail_details) < 25:
                    fail_details.append(
                        f"  FAIL [EQ] {loc_name:12s} JD={jd:.1f}: "
                        f'dRA={dra:.2f}" dDec={ddec:.2f}"'
                    )

    print(
        f"  Part 2: {p2_pass}/{p2_pass + p2_fail + p2_err} passed, {p2_fail} fail, {p2_err} err"
    )

    # =========================================================================
    # PART 3: Topocentric vs geocentric parallax magnitude check
    # =========================================================================
    print("\n--- PART 3: Parallax Magnitude Verification ---")
    p3_pass = 0
    p3_fail = 0
    p3_err = 0

    for loc_lat, loc_lon, loc_alt, loc_name in LOCATIONS[:8]:
        swe.set_topo(loc_lon, loc_lat, loc_alt)
        ephem.set_topo(loc_lon, loc_lat, loc_alt)

        for jd in EPOCHS[:6]:
            total += 1
            flags_geo = FLG_SWIEPH | FLG_SPEED
            flags_topo = FLG_SWIEPH | FLG_SPEED | FLG_TOPOCTR

            try:
                le_geo = ephem.calc_ut(jd, MOON, flags_geo)[0]
                le_topo = ephem.calc_ut(jd, MOON, flags_topo)[0]
                se_geo = swe.calc_ut(jd, MOON, flags_geo)[0]
                se_topo = swe.calc_ut(jd, MOON, flags_topo)[0]
            except Exception as e:
                p3_err += 1
                errors += 1
                continue

            # Parallax in longitude
            le_parallax_lon = abs(le_topo[0] - le_geo[0]) * 3600
            se_parallax_lon = abs(se_topo[0] - se_geo[0]) * 3600

            # The Moon's horizontal parallax is ~3400" (57'), so topocentric
            # correction in longitude should be 0-3600" depending on geometry
            parallax_diff = abs(le_parallax_lon - se_parallax_lon)

            # Also check that parallax is in reasonable range
            max_parallax = 4000  # arcsec (generous bound)

            if parallax_diff <= 3.0 and le_parallax_lon <= max_parallax:
                p3_pass += 1
                passed += 1
            else:
                p3_fail += 1
                failed += 1
                if len(fail_details) < 30:
                    fail_details.append(
                        f"  FAIL [PAR] {loc_name:12s} JD={jd:.1f}: "
                        f'LE_par={le_parallax_lon:.1f}" SE_par={se_parallax_lon:.1f}" '
                        f'diff={parallax_diff:.1f}"'
                    )

    print(
        f"  Part 3: {p3_pass}/{p3_pass + p3_fail + p3_err} passed, {p3_fail} fail, {p3_err} err"
    )

    # =========================================================================
    # PART 4: Topocentric Sun (should have very small parallax ~8.8")
    # =========================================================================
    print("\n--- PART 4: Topocentric Sun (small parallax check) ---")
    p4_pass = 0
    p4_fail = 0
    p4_err = 0

    SUN = 0

    for loc_lat, loc_lon, loc_alt, loc_name in LOCATIONS[:8]:
        swe.set_topo(loc_lon, loc_lat, loc_alt)
        ephem.set_topo(loc_lon, loc_lat, loc_alt)

        for jd in EPOCHS[:6]:
            total += 1
            flags = FLG_SWIEPH | FLG_SPEED | FLG_TOPOCTR

            try:
                se_pos = swe.calc_ut(jd, SUN, flags)[0]
                le_pos = ephem.calc_ut(jd, SUN, flags)[0]
            except Exception as e:
                p4_err += 1
                errors += 1
                continue

            dlon = abs(se_pos[0] - le_pos[0]) * 3600
            dlat = abs(se_pos[1] - le_pos[1]) * 3600

            if dlon > 180 * 3600:
                dlon = 360 * 3600 - dlon

            if dlon <= 1.0 and dlat <= 1.0:
                p4_pass += 1
                passed += 1
            else:
                p4_fail += 1
                failed += 1
                if len(fail_details) < 35:
                    fail_details.append(
                        f"  FAIL [SUN] {loc_name:12s} JD={jd:.1f}: "
                        f'dLon={dlon:.3f}" dLat={dlat:.3f}"'
                    )

    print(
        f"  Part 4: {p4_pass}/{p4_pass + p4_fail + p4_err} passed, {p4_fail} fail, {p4_err} err"
    )

    # =========================================================================
    # PART 5: Topocentric planets (Mars, Jupiter, Saturn)
    # =========================================================================
    print("\n--- PART 5: Topocentric Planets (Mars/Jupiter/Saturn) ---")
    p5_pass = 0
    p5_fail = 0
    p5_err = 0

    PLANETS = [(4, "Mars"), (5, "Jupiter"), (6, "Saturn")]

    for loc_lat, loc_lon, loc_alt, loc_name in LOCATIONS[:5]:
        swe.set_topo(loc_lon, loc_lat, loc_alt)
        ephem.set_topo(loc_lon, loc_lat, loc_alt)

        for body_id, body_name in PLANETS:
            for jd in EPOCHS[:6]:
                total += 1
                flags = FLG_SWIEPH | FLG_SPEED | FLG_TOPOCTR

                try:
                    se_pos = swe.calc_ut(jd, body_id, flags)[0]
                    le_pos = ephem.calc_ut(jd, body_id, flags)[0]
                except Exception as e:
                    p5_err += 1
                    errors += 1
                    continue

                dlon = abs(se_pos[0] - le_pos[0]) * 3600
                dlat = abs(se_pos[1] - le_pos[1]) * 3600

                if dlon > 180 * 3600:
                    dlon = 360 * 3600 - dlon

                if dlon <= 1.0 and dlat <= 1.0:
                    p5_pass += 1
                    passed += 1
                else:
                    p5_fail += 1
                    failed += 1
                    if len(fail_details) < 40:
                        fail_details.append(
                            f"  FAIL [{body_name}] {loc_name:12s} JD={jd:.1f}: "
                            f'dLon={dlon:.3f}" dLat={dlat:.3f}"'
                        )

    print(
        f"  Part 5: {p5_pass}/{p5_pass + p5_fail + p5_err} passed, {p5_fail} fail, {p5_err} err"
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
        f"ROUND 97 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
