#!/usr/bin/env python3
"""
Round 9: Deep Topocentric Position Audit
==========================================
Compares libephemeris topocentric calculations against pyswisseph:
  P1: Topocentric planet positions — all major bodies at multiple locations
  P2: Topocentric Moon parallax — largest parallactic shift (~1°)
  P3: Topocentric + equatorial combined flags
  P4: Topocentric + sidereal combined flags
  P5: Topocentric + J2000 combined flags
  P6: Topocentric at extreme altitudes (sea level, mountain, airplane)
  P7: Topocentric velocity consistency
  P8: Topocentric Sun position (solar parallax ~8.8")
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


# ============================================================================
# CONSTANTS
# ============================================================================

FLG_SPEED = swe.FLG_SPEED
FLG_TOPOCTR = swe.FLG_TOPOCTR
FLG_EQUATORIAL = swe.FLG_EQUATORIAL
FLG_SIDEREAL = swe.FLG_SIDEREAL
FLG_J2000 = 32  # FLG_J2000

# Observer locations: (name, lon, lat, alt_m)
LOCATIONS = [
    ("Rome", 12.4964, 41.9028, 21.0),
    ("NewYork", -74.0060, 40.7128, 10.0),
    ("Tokyo", 139.6917, 35.6895, 40.0),
    ("Sydney", 151.2093, -33.8688, 58.0),
    ("NorthPole", 0.0, 89.0, 0.0),
    ("Quito", -78.4678, -0.1807, 2850.0),  # High altitude equator
]

# Bodies to test
BODIES = [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MERCURY, "Mercury"),
    (swe.VENUS, "Venus"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
    (swe.URANUS, "Uranus"),
    (swe.NEPTUNE, "Neptune"),
    (swe.PLUTO, "Pluto"),
    (swe.MEAN_NODE, "MeanNode"),
    (swe.TRUE_NODE, "TrueNode"),
    (swe.MEAN_APOG, "MeanLilith"),
]

# Test epochs
JD_2024 = 2460310.5  # 2024-Jan-15
JD_2000 = 2451545.0  # J2000.0
JD_1990 = 2448000.5  # 1990


# ============================================================================
# PART 1: Topocentric planet positions at multiple locations
# ============================================================================


def test_part1_topo_positions():
    print("\n" + "=" * 70)
    print("PART 1: Topocentric Planet Positions — Multiple Locations")
    print("=" * 70)

    jd = JD_2024
    flags = FLG_SPEED | FLG_TOPOCTR

    for loc_name, lon, lat, alt in LOCATIONS:
        swe.set_topo(lon, lat, alt)
        ephem.set_topo(lon, lat, alt)

        for body_id, body_name in BODIES:
            test_name = f"P1/{loc_name}/{body_name}"
            try:
                pos_se, retflag_se = swe.calc_ut(jd, body_id, flags)
                pos_le, retflag_le = ephem.calc_ut(jd, body_id, flags)

                lon_se = float(pos_se[0])
                lon_le = float(pos_le[0])
                lat_se = float(pos_se[1])
                lat_le = float(pos_le[1])
                speed_se = float(pos_se[3])
                speed_le = float(pos_le[3])

                lon_diff = abs(lon_se - lon_le)
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff
                lat_diff = abs(lat_se - lat_le)
                speed_diff = abs(speed_se - speed_le)

                # Topocentric: allow up to 1" for planets, 2" for Moon (parallax)
                lon_tol = 0.001 if body_name == "Moon" else 0.0005  # degrees
                lat_tol = 0.001 if body_name == "Moon" else 0.0005
                speed_tol = 0.001  # deg/day

                ok_lon = lon_diff < lon_tol
                ok_lat = lat_diff < lat_tol
                ok_speed = speed_diff < speed_tol

                record(
                    f"{test_name}/lon",
                    ok_lon,
                    f"SE={lon_se:.6f} LE={lon_le:.6f} "
                    f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
                )
                record(
                    f"{test_name}/lat",
                    ok_lat,
                    f"SE={lat_se:.6f} LE={lat_le:.6f} diff={lat_diff:.6f}°",
                )
                record(
                    f"{test_name}/speed",
                    ok_speed,
                    f"SE={speed_se:.6f} LE={speed_le:.6f} diff={speed_diff:.6f}°/day",
                )

            except Exception as e:
                record(test_name, False, f"ERROR: {e}")
                traceback.print_exc()


# ============================================================================
# PART 2: Moon parallax verification — topocentric shift should be ~0.5-1°
# ============================================================================


def test_part2_moon_parallax():
    print("\n" + "=" * 70)
    print("PART 2: Moon Parallax — Topocentric vs Geocentric Shift")
    print("=" * 70)

    jd = JD_2024

    for loc_name, lon, lat, alt in LOCATIONS[:3]:  # Test 3 locations
        swe.set_topo(lon, lat, alt)
        ephem.set_topo(lon, lat, alt)

        test_name = f"P2/{loc_name}/Moon"
        try:
            # Geocentric
            geo_se, _ = swe.calc_ut(jd, swe.MOON, FLG_SPEED)
            geo_le, _ = ephem.calc_ut(jd, swe.MOON, FLG_SPEED)

            # Topocentric
            topo_se, _ = swe.calc_ut(jd, swe.MOON, FLG_SPEED | FLG_TOPOCTR)
            topo_le, _ = ephem.calc_ut(jd, swe.MOON, FLG_SPEED | FLG_TOPOCTR)

            # Parallax shift
            shift_se = abs(float(geo_se[0]) - float(topo_se[0]))
            if shift_se > 180:
                shift_se = 360 - shift_se
            shift_le = abs(float(geo_le[0]) - float(topo_le[0]))
            if shift_le > 180:
                shift_le = 360 - shift_le

            shift_diff = abs(shift_se - shift_le)

            # Moon parallax should be 0.5-1.0° — shifts should agree within 0.001°
            ok = shift_diff < 0.001
            record(
                f"{test_name}/parallax",
                ok,
                f"SE_shift={shift_se:.6f}° LE_shift={shift_le:.6f}° "
                f'diff={shift_diff:.6f}° ({shift_diff * 3600:.2f}")',
            )

            # Verify parallax is in reasonable range (0.3° - 1.1°)
            ok_range = 0.3 < shift_le < 1.1
            record(
                f"{test_name}/range",
                ok_range,
                f"LE_shift={shift_le:.6f}° (expect 0.3°-1.1°)",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 3: Topocentric + Equatorial combined flags
# ============================================================================


def test_part3_topo_equatorial():
    print("\n" + "=" * 70)
    print("PART 3: Topocentric + Equatorial Combined Flags")
    print("=" * 70)

    jd = JD_2024
    flags = FLG_SPEED | FLG_TOPOCTR | FLG_EQUATORIAL

    loc_name, lon, lat, alt = LOCATIONS[0]  # Rome
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)

    bodies_subset = BODIES[:7]  # Sun through Saturn

    for body_id, body_name in bodies_subset:
        test_name = f"P3/topo+eq/{body_name}"
        try:
            pos_se, _ = swe.calc_ut(jd, body_id, flags)
            pos_le, _ = ephem.calc_ut(jd, body_id, flags)

            ra_se = float(pos_se[0])
            ra_le = float(pos_le[0])
            dec_se = float(pos_se[1])
            dec_le = float(pos_le[1])

            ra_diff = abs(ra_se - ra_le)
            if ra_diff > 180:
                ra_diff = 360 - ra_diff
            dec_diff = abs(dec_se - dec_le)

            tol = 0.001 if body_name == "Moon" else 0.0005

            ok_ra = ra_diff < tol
            ok_dec = dec_diff < tol

            record(
                f"{test_name}/RA",
                ok_ra,
                f"SE={ra_se:.6f} LE={ra_le:.6f} "
                f'diff={ra_diff:.6f}° ({ra_diff * 3600:.2f}")',
            )
            record(
                f"{test_name}/Dec",
                ok_dec,
                f"SE={dec_se:.6f} LE={dec_le:.6f} diff={dec_diff:.6f}°",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 4: Topocentric + Sidereal combined flags
# ============================================================================


def test_part4_topo_sidereal():
    print("\n" + "=" * 70)
    print("PART 4: Topocentric + Sidereal Combined Flags")
    print("=" * 70)

    jd = JD_2024
    swe.set_sid_mode(1)  # Lahiri
    ephem.set_sid_mode(1)
    flags = FLG_SPEED | FLG_TOPOCTR | FLG_SIDEREAL

    loc_name, lon, lat, alt = LOCATIONS[0]  # Rome
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)

    bodies_subset = BODIES[:7]

    for body_id, body_name in bodies_subset:
        test_name = f"P4/topo+sid/{body_name}"
        try:
            pos_se, _ = swe.calc_ut(jd, body_id, flags)
            pos_le, _ = ephem.calc_ut(jd, body_id, flags)

            lon_se = float(pos_se[0])
            lon_le = float(pos_le[0])
            lon_diff = abs(lon_se - lon_le)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            tol = 0.001 if body_name == "Moon" else 0.0005

            ok = lon_diff < tol
            record(
                f"{test_name}/lon",
                ok,
                f"SE={lon_se:.6f} LE={lon_le:.6f} "
                f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0)


# ============================================================================
# PART 5: Topocentric + J2000 combined flags
# ============================================================================


def test_part5_topo_j2000():
    print("\n" + "=" * 70)
    print("PART 5: Topocentric + J2000 Combined Flags")
    print("=" * 70)

    jd = JD_2024
    flags = FLG_SPEED | FLG_TOPOCTR | FLG_J2000

    loc_name, lon, lat, alt = LOCATIONS[0]  # Rome
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)

    bodies_subset = BODIES[:7]

    for body_id, body_name in bodies_subset:
        test_name = f"P5/topo+J2000/{body_name}"
        try:
            pos_se, _ = swe.calc_ut(jd, body_id, flags)
            pos_le, _ = ephem.calc_ut(jd, body_id, flags)

            lon_se = float(pos_se[0])
            lon_le = float(pos_le[0])
            lon_diff = abs(lon_se - lon_le)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            tol = 0.001 if body_name == "Moon" else 0.0005

            ok = lon_diff < tol
            record(
                f"{test_name}/lon",
                ok,
                f"SE={lon_se:.6f} LE={lon_le:.6f} "
                f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 6: Topocentric at extreme altitudes
# ============================================================================


def test_part6_extreme_altitudes():
    print("\n" + "=" * 70)
    print("PART 6: Topocentric at Extreme Altitudes")
    print("=" * 70)

    jd = JD_2024
    flags = FLG_SPEED | FLG_TOPOCTR
    lon, lat = 12.4964, 41.9028  # Rome coordinates

    altitudes = [
        ("SeaLevel", 0.0),
        ("City", 100.0),
        ("Mountain", 3500.0),
        ("HighMountain", 5500.0),
        ("Airplane", 11000.0),
    ]

    for alt_name, alt in altitudes:
        swe.set_topo(lon, lat, alt)
        ephem.set_topo(lon, lat, alt)

        # Test Moon (most affected by altitude)
        test_name = f"P6/{alt_name}/Moon"
        try:
            pos_se, _ = swe.calc_ut(jd, swe.MOON, flags)
            pos_le, _ = ephem.calc_ut(jd, swe.MOON, flags)

            lon_se = float(pos_se[0])
            lon_le = float(pos_le[0])
            lon_diff = abs(lon_se - lon_le)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            ok = lon_diff < 0.001
            record(
                f"{test_name}/lon",
                ok,
                f"SE={lon_se:.6f} LE={lon_le:.6f} "
                f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")

        # Test Sun
        test_name = f"P6/{alt_name}/Sun"
        try:
            pos_se, _ = swe.calc_ut(jd, swe.SUN, flags)
            pos_le, _ = ephem.calc_ut(jd, swe.SUN, flags)

            lon_se = float(pos_se[0])
            lon_le = float(pos_le[0])
            lon_diff = abs(lon_se - lon_le)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            ok = lon_diff < 0.0005
            record(
                f"{test_name}/lon",
                ok,
                f"SE={lon_se:.6f} LE={lon_le:.6f} "
                f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 7: Topocentric velocity consistency
# ============================================================================


def test_part7_topo_velocity():
    print("\n" + "=" * 70)
    print("PART 7: Topocentric Velocity Consistency")
    print("=" * 70)

    jd = JD_2024
    flags = FLG_SPEED | FLG_TOPOCTR

    loc_name, lon, lat, alt = LOCATIONS[0]  # Rome
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)

    bodies_subset = BODIES[:7]  # Sun through Saturn

    for body_id, body_name in bodies_subset:
        test_name = f"P7/topo_vel/{body_name}"
        try:
            pos_se, _ = swe.calc_ut(jd, body_id, flags)
            pos_le, _ = ephem.calc_ut(jd, body_id, flags)

            speed_lon_se = float(pos_se[3])
            speed_lon_le = float(pos_le[3])
            speed_lat_se = float(pos_se[4])
            speed_lat_le = float(pos_le[4])

            speed_lon_diff = abs(speed_lon_se - speed_lon_le)
            speed_lat_diff = abs(speed_lat_se - speed_lat_le)

            # Topocentric speeds can differ more than geocentric due to
            # diurnal motion parallax — allow 0.005°/day for Moon, 0.001 others
            lon_tol = 0.005 if body_name == "Moon" else 0.001
            lat_tol = 0.005 if body_name == "Moon" else 0.001

            ok_lon = speed_lon_diff < lon_tol
            ok_lat = speed_lat_diff < lat_tol

            record(
                f"{test_name}/speed_lon",
                ok_lon,
                f"SE={speed_lon_se:.6f} LE={speed_lon_le:.6f} "
                f"diff={speed_lon_diff:.6f}°/day",
            )
            record(
                f"{test_name}/speed_lat",
                ok_lat,
                f"SE={speed_lat_se:.6f} LE={speed_lat_le:.6f} "
                f"diff={speed_lat_diff:.6f}°/day",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# PART 8: Topocentric Sun — solar parallax ~8.8"
# ============================================================================


def test_part8_topo_sun():
    print("\n" + "=" * 70)
    print("PART 8: Topocentric Sun Position at Multiple Epochs")
    print("=" * 70)

    flags = FLG_SPEED | FLG_TOPOCTR
    epochs = [
        ("J2000", JD_2000),
        ("2024-Jan", JD_2024),
        ("1990", JD_1990),
    ]

    loc_name, lon, lat, alt = LOCATIONS[0]  # Rome
    swe.set_topo(lon, lat, alt)
    ephem.set_topo(lon, lat, alt)

    for epoch_name, jd in epochs:
        test_name = f"P8/Sun/{epoch_name}"
        try:
            pos_se, _ = swe.calc_ut(jd, swe.SUN, flags)
            pos_le, _ = ephem.calc_ut(jd, swe.SUN, flags)

            lon_se = float(pos_se[0])
            lon_le = float(pos_le[0])
            lat_se = float(pos_se[1])
            lat_le = float(pos_le[1])
            dist_se = float(pos_se[2])
            dist_le = float(pos_le[2])

            lon_diff = abs(lon_se - lon_le)
            if lon_diff > 180:
                lon_diff = 360 - lon_diff
            lat_diff = abs(lat_se - lat_le)
            dist_diff = abs(dist_se - dist_le)

            ok_lon = lon_diff < 0.0005
            ok_lat = lat_diff < 0.0005
            ok_dist = dist_diff < 0.0001  # AU

            record(
                f"{test_name}/lon",
                ok_lon,
                f"SE={lon_se:.6f} LE={lon_le:.6f} "
                f'diff={lon_diff:.6f}° ({lon_diff * 3600:.2f}")',
            )
            record(
                f"{test_name}/lat",
                ok_lat,
                f"SE={lat_se:.6f} LE={lat_le:.6f} diff={lat_diff:.6f}°",
            )
            record(
                f"{test_name}/dist",
                ok_dist,
                f"SE={dist_se:.8f} LE={dist_le:.8f} diff={dist_diff:.8f} AU",
            )

        except Exception as e:
            record(test_name, False, f"ERROR: {e}")


# ============================================================================
# MAIN
# ============================================================================


def main():
    global total, passed, failed, skipped, failures

    print("=" * 70)
    print("ROUND 9: Deep Topocentric Position Audit")
    print("=" * 70)
    t0 = time.time()

    test_part1_topo_positions()
    test_part2_moon_parallax()
    test_part3_topo_equatorial()
    test_part4_topo_sidereal()
    test_part5_topo_j2000()
    test_part6_extreme_altitudes()
    test_part7_topo_velocity()
    test_part8_topo_sun()

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

    # Reset topo
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
