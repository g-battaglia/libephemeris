#!/usr/bin/env python3
"""Round 79: Crossing Functions Deep Sweep

Tests all crossing functions comprehensively:
- solcross_ut: Sun crossing specific longitudes
- mooncross_ut: Moon crossing specific longitudes
- mooncross_node_ut: Moon crossing its own nodes
- cross_ut: Generic planet crossing longitudes
- helio_cross_ut: Heliocentric crossing
- find_station_ut: Planetary stations (retrograde/direct)
"""

from __future__ import annotations

import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0


def check(label, se_val, le_val, tol_sec, is_time=False):
    global passed, failed
    if is_time:
        # Time comparison in seconds
        diff = abs(se_val - le_val) * 86400.0  # days to seconds
        if diff < tol_sec:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL {label}: SE={se_val:.8f} LE={le_val:.8f} diff={diff:.2f}s")
    else:
        diff = abs(se_val - le_val) * 3600.0  # degrees to arcsec
        if diff < tol_sec:
            passed += 1
        else:
            failed += 1
            print(f'  FAIL {label}: SE={se_val:.6f} LE={le_val:.6f} diff={diff:.2f}"')


print("=" * 70)
print("ROUND 79: Crossing Functions Deep Sweep")
print("=" * 70)

# ============================================================
# P1: Sun crossing specific longitudes (ingresses)
# ============================================================
print("\n=== P1: Sun crossing longitudes (solcross_ut) ===")

# Test Sun crossing every 30° (zodiac ingresses) over 10 years
jd_start = 2451545.0  # J2000.0
for year_offset in range(10):
    jd = jd_start + year_offset * 365.25
    for target_lon in [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]:
        label = f"solcross y+{year_offset} lon={target_lon}"
        try:
            se_jd = swe.solcross_ut(target_lon, jd, swe.FLG_SWIEPH)
            le_jd = ephem.solcross_ut(target_lon, jd, 2)
            # Verify the Sun is actually at the target longitude
            check(label, se_jd, le_jd, 2.0, is_time=True)  # 2 second tolerance
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

p1_p, p1_f = passed, failed
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Moon crossing specific longitudes
# ============================================================
print("\n=== P2: Moon crossing longitudes (mooncross_ut) ===")

# Moon crosses each longitude ~once per month
for month in range(24):  # 2 years monthly
    jd = jd_start + month * 29.53
    for target_lon in [0, 45, 90, 135, 180, 225, 270, 315]:
        label = f"mooncross m={month} lon={target_lon}"
        try:
            se_jd = swe.mooncross_ut(target_lon, jd, swe.FLG_SWIEPH)
            le_jd = ephem.mooncross_ut(target_lon, jd, 2)
            check(label, se_jd, le_jd, 1.0, is_time=True)  # 1 second tolerance
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

p2_p, p2_f = passed - p1_p, failed - p1_f
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Moon crossing its own nodes
# ============================================================
print("\n=== P3: Moon node crossings (mooncross_node_ut) ===")

# Moon crosses nodes every ~13.6 days
for i in range(50):  # ~50 node crossings over ~2 years
    jd = jd_start + i * 13.6
    label = f"mooncross_node i={i}"
    try:
        se_result = swe.mooncross_node_ut(jd, swe.FLG_SWIEPH)
        # se_result is (jd_cross, xlon, xlat, xdist) or similar
        se_jd = se_result[0]
        se_lon = se_result[1]
        se_lat = se_result[2]

        le_result = ephem.mooncross_node_ut(jd, 2)
        le_jd = le_result[0]
        le_lon = le_result[1]
        le_lat = le_result[2]

        # Known: SE returns TT as UT (confirmed in Round 4)
        # so we expect ~69s systematic offset
        diff_s = abs(se_jd - le_jd) * 86400.0
        if diff_s < 120.0:  # Allow up to 120s (known ~69s TT/UT offset)
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL {label} time: SE={se_jd:.8f} LE={le_jd:.8f} diff={diff_s:.1f}s"
            )

        # Check that latitude at crossing is near zero
        if abs(le_lat) < 0.01:  # latitude should be very small at node
            passed += 1
        else:
            failed += 1
            print(f"  FAIL {label} lat: {le_lat:.6f}° (should be ~0)")

    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

p3_p, p3_f = passed - p1_p - p2_p, failed - p1_f - p2_f
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Generic planet crossing (cross_ut)
# ============================================================
print("\n=== P4: Planet crossing longitudes (cross_ut) ===")

planets = [
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
]

for body_id, name in planets:
    for year_offset in range(5):
        jd = jd_start + year_offset * 365.25
        # Get current planet position to pick a target ahead
        se_pos = swe.calc_ut(jd, body_id, swe.FLG_SWIEPH | swe.FLG_SPEED)
        current_lon = se_pos[0][0]
        # Target: 30° ahead
        target_lon = (current_lon + 30.0) % 360.0
        label = f"cross_ut {name} y+{year_offset} lon={target_lon:.1f}"
        try:
            se_jd = swe.cross_ut(body_id, target_lon, jd, swe.FLG_SWIEPH)
            le_jd = ephem.cross_ut(body_id, target_lon, jd, 2)
            check(label, se_jd, le_jd, 5.0, is_time=True)  # 5 second tolerance
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Heliocentric crossing
# ============================================================
print("\n=== P5: Heliocentric crossing (helio_cross_ut) ===")

for body_id, name in planets:
    for year_offset in range(5):
        jd = jd_start + year_offset * 365.25
        # Get heliocentric position
        se_pos = swe.calc_ut(jd, body_id, swe.FLG_SWIEPH | swe.FLG_HELCTR)
        current_lon = se_pos[0][0]
        target_lon = (current_lon + 30.0) % 360.0
        label = f"helio_cross {name} y+{year_offset} lon={target_lon:.1f}"
        try:
            se_jd = swe.helio_cross_ut(body_id, target_lon, jd, swe.FLG_SWIEPH)
            le_jd = ephem.helio_cross_ut(body_id, target_lon, jd, 2)
            check(
                label, se_jd, le_jd, 10.0, is_time=True
            )  # 10 second tolerance (helio slower convergence)
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: Planetary stations (find_station_ut) — direct & retrograde
# ============================================================
print("\n=== P6: Planetary stations (find_station_ut) ===")

station_planets = [
    (swe.MERCURY, "Mercury"),
    (swe.VENUS, "Venus"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
]

for body_id, name in station_planets:
    for year_offset in range(5):
        jd = jd_start + year_offset * 365.25
        for direction in [0, 1]:  # 0=next retrograde station, 1=next direct station
            dir_name = "retro" if direction == 0 else "direct"
            label = f"station {name} y+{year_offset} {dir_name}"
            try:
                # SE: find_station_ut(body, jd, flags, backward)
                # backward: 0=forward, 1=backward
                # direction: we use backward=0 always, just look forward
                se_result = None
                try:
                    se_result = swe.find_station_ut(
                        body_id, jd, swe.FLG_SWIEPH, direction
                    )
                except Exception:
                    pass

                if se_result is None or se_result == 0:
                    continue

                # LE API: find_station_ut(body, jd_start, flags, direction)
                # direction: 0=retrograde station, 1=direct station
                le_result = ephem.find_station_ut(body_id, jd, 2, direction)

                if isinstance(se_result, (int, float)):
                    se_jd = float(se_result)
                elif isinstance(se_result, tuple):
                    se_jd = se_result[0]
                else:
                    continue

                if isinstance(le_result, (int, float)):
                    le_jd = float(le_result)
                elif isinstance(le_result, tuple):
                    le_jd = le_result[0]
                else:
                    continue

                check(
                    label, se_jd, le_jd, 60.0, is_time=True
                )  # 60s tolerance for stations

            except Exception as e:
                errors += 1
                print(f"  ERROR {label}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Sun crossing verification — check actual position at crossing
# ============================================================
print("\n=== P7: Position at crossing verification ===")

for target_lon in [0, 90, 180, 270]:
    for year_offset in range(10):
        jd = jd_start + year_offset * 365.25
        label = f"solcross_verify y+{year_offset} lon={target_lon}"
        try:
            le_jd = ephem.solcross_ut(target_lon, jd, 2)
            # Verify Sun is actually at target longitude
            le_pos = ephem.calc_ut(le_jd, 0)  # Sun
            actual_lon = le_pos[0]
            diff_arcsec = abs(actual_lon - target_lon) * 3600.0
            # Handle wraparound
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec
            if diff_arcsec < 0.001:  # sub-milliarcsecond precision
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: target={target_lon} actual={actual_lon:.8f} diff={diff_arcsec:.6f}"'
                )
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P8: Moon crossing verification — check actual position at crossing
# ============================================================
print("\n=== P8: Moon position at crossing verification ===")

for month in range(12):
    jd = jd_start + month * 29.53
    for target_lon in [0, 90, 180, 270]:
        label = f"mooncross_verify m={month} lon={target_lon}"
        try:
            le_jd = ephem.mooncross_ut(target_lon, jd, 2)
            le_pos = ephem.calc_ut(le_jd, 1)  # Moon
            actual_lon = le_pos[0]
            diff_arcsec = abs(actual_lon - target_lon) * 3600.0
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec
            if diff_arcsec < 0.01:  # 10 mas precision
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: target={target_lon} actual={actual_lon:.8f} diff={diff_arcsec:.6f}"'
                )
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P9: Cross_ut with inner planets (fast-moving)
# ============================================================
print("\n=== P9: Inner planet crossings (Mercury, Venus) ===")

inner_planets = [
    (swe.MERCURY, "Mercury"),
    (swe.VENUS, "Venus"),
]

for body_id, name in inner_planets:
    for month in range(12):
        jd = jd_start + month * 30.0
        se_pos = swe.calc_ut(jd, body_id, swe.FLG_SWIEPH | swe.FLG_SPEED)
        current_lon = se_pos[0][0]
        # Skip if retrograde (crossing may be ambiguous)
        if se_pos[0][3] < 0:
            continue
        target_lon = (current_lon + 20.0) % 360.0
        label = f"cross_ut {name} m={month} lon={target_lon:.1f}"
        try:
            se_jd = swe.cross_ut(body_id, target_lon, jd, swe.FLG_SWIEPH)
            le_jd = ephem.cross_ut(body_id, target_lon, jd, 2)
            check(label, se_jd, le_jd, 5.0, is_time=True)
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P10: Edge case — crossing at 0°/360° boundary
# ============================================================
print("\n=== P10: 0°/360° boundary crossings ===")

for year_offset in range(5):
    jd = jd_start + year_offset * 365.25
    # Sun crosses 0° (Aries ingress)
    label = f"solcross_0deg y+{year_offset}"
    try:
        se_jd = swe.solcross_ut(0.0, jd, swe.FLG_SWIEPH)
        le_jd = ephem.solcross_ut(0.0, jd, 2)
        check(label, se_jd, le_jd, 2.0, is_time=True)
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

    # Moon crossing 359.99° (near boundary)
    label = f"mooncross_359.99 y+{year_offset}"
    try:
        se_jd = swe.mooncross_ut(359.99, jd, swe.FLG_SWIEPH)
        le_jd = ephem.mooncross_ut(359.99, jd, 2)
        check(label, se_jd, le_jd, 1.0, is_time=True)
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

print(f"  After P10: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# Summary
# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 79 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
