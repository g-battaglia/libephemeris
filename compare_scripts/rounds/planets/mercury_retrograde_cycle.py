#!/usr/bin/env python3
"""Round 56: Mercury Retrograde Precision

Verify Mercury retrograde periods:
- Mercury retrogrades ~3 times per year, lasting ~24 days each
- Station (retrograde) occurs when lon_speed = 0 and transitions negative
- Compare station timing, position at station, and speed reversal precision
- Also verify Mars retrograde (~every 26 months) for slower planet comparison
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

MERCURY = 2
MARS = 4
JUPITER = 5
SATURN = 6
FLAGS = 256  # FLG_SPEED

print("=" * 70)
print("ROUND 56: Mercury Retrograde Precision")
print("=" * 70)


def find_stations(body_id, jd_start, jd_end, step=1.0):
    """Find station points (speed = 0 crossing) for a body using LE."""
    stations = []
    jd = jd_start
    prev_speed = None
    while jd < jd_end:
        try:
            result = ephem.calc_ut(jd, body_id, FLAGS)
            speed = result[0][3]  # lon speed deg/day
            if prev_speed is not None:
                if prev_speed > 0 and speed < 0:
                    # Station retrograde - bisect to find exact time
                    stations.append(("R", _bisect_station(body_id, jd - step, jd)))
                elif prev_speed < 0 and speed > 0:
                    # Station direct
                    stations.append(("D", _bisect_station(body_id, jd - step, jd)))
            prev_speed = speed
        except Exception:
            pass
        jd += step
    return stations


def _bisect_station(body_id, jd_lo, jd_hi, tol=1e-8):
    """Bisect to find when lon speed = 0."""
    for _ in range(60):
        jd_mid = (jd_lo + jd_hi) / 2
        result = ephem.calc_ut(jd_mid, body_id, FLAGS)
        speed = result[0][3]
        if abs(speed) < 1e-10:
            return jd_mid
        # Get speed at lo
        lo_result = ephem.calc_ut(jd_lo, body_id, FLAGS)
        lo_speed = lo_result[0][3]
        if (lo_speed > 0 and speed > 0) or (lo_speed < 0 and speed < 0):
            jd_lo = jd_mid
        else:
            jd_hi = jd_mid
        if jd_hi - jd_lo < tol:
            break
    return (jd_lo + jd_hi) / 2


def find_stations_se(body_id, jd_start, jd_end, step=1.0):
    """Find station points using SE."""
    stations = []
    jd = jd_start
    prev_speed = None
    while jd < jd_end:
        try:
            result = swe.calc_ut(jd, body_id, FLAGS)
            speed = result[0][3]
            if prev_speed is not None:
                if prev_speed > 0 and speed < 0:
                    stations.append(("R", _bisect_station_se(body_id, jd - step, jd)))
                elif prev_speed < 0 and speed > 0:
                    stations.append(("D", _bisect_station_se(body_id, jd - step, jd)))
            prev_speed = speed
        except Exception:
            pass
        jd += step
    return stations


def _bisect_station_se(body_id, jd_lo, jd_hi, tol=1e-8):
    """Bisect to find when lon speed = 0 using SE."""
    for _ in range(60):
        jd_mid = (jd_lo + jd_hi) / 2
        result = swe.calc_ut(jd_mid, body_id, FLAGS)
        speed = result[0][3]
        if abs(speed) < 1e-10:
            return jd_mid
        lo_result = swe.calc_ut(jd_lo, body_id, FLAGS)
        lo_speed = lo_result[0][3]
        if (lo_speed > 0 and speed > 0) or (lo_speed < 0 and speed < 0):
            jd_lo = jd_mid
        else:
            jd_hi = jd_mid
        if jd_hi - jd_lo < tol:
            break
    return (jd_lo + jd_hi) / 2


# ============================================================
# P1: Mercury stations 2020-2030 (timing comparison)
# ============================================================
print("\n=== P1: Mercury station timing 2020-2030 ===")

jd_2020 = swe.julday(2020, 1, 1, 0.0)
jd_2030 = swe.julday(2030, 1, 1, 0.0)

le_stations = find_stations(MERCURY, jd_2020, jd_2030, step=0.5)
se_stations = find_stations_se(MERCURY, jd_2020, jd_2030, step=0.5)

print(
    f"  LE found {len(le_stations)} Mercury stations, SE found {len(se_stations)} stations"
)

# Match stations and compare timing
if len(le_stations) == len(se_stations):
    passed += 1
    for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(
        zip(le_stations, se_stations)
    ):
        if le_type == se_type:
            diff_sec = abs(le_jd - se_jd) * 86400  # seconds
            # Station timing: should agree within 60 seconds
            if diff_sec < 60.0:
                passed += 1
            else:
                failed += 1
                # Convert to date for readability
                y, m, d, h = swe.revjul(le_jd)
                print(
                    f"  FAIL P1 station {i} ({le_type}) {y}-{m:02d}-{d:02d}: diff={diff_sec:.1f}s"
                )
        else:
            failed += 1
            print(f"  FAIL P1 station {i}: LE type={le_type} SE type={se_type}")
else:
    failed += 1
    print(
        f"  FAIL P1: station count mismatch LE={len(le_stations)} SE={len(se_stations)}"
    )

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Mercury position at station points
# ============================================================
print("\n=== P2: Mercury longitude at station points ===")

for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(zip(le_stations, se_stations)):
    try:
        le_result = ephem.calc_ut(le_jd, MERCURY, FLAGS)
        se_result = swe.calc_ut(se_jd, MERCURY, FLAGS)

        le_lon = le_result[0][0]
        se_lon = se_result[0][0]

        diff_arcsec = abs(le_lon - se_lon) * 3600
        if diff_arcsec > 180 * 3600:
            diff_arcsec = 360 * 3600 - diff_arcsec

        # Position at station should agree within 1"
        if diff_arcsec < 1.0:
            passed += 1
        else:
            y, m, d, h = swe.revjul(le_jd)
            failed += 1
            print(
                f'  FAIL P2 station {i} ({le_type}) {y}-{m:02d}: SE={se_lon:.6f}° LE={le_lon:.6f}° diff={diff_arcsec:.2f}"'
            )

    except Exception as e:
        errors += 1
        print(f"  ERROR P2 station {i}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Mercury retrograde duration (R to D station interval)
# ============================================================
print("\n=== P3: Mercury retrograde duration ===")

for i in range(0, len(le_stations) - 1, 2):
    if i + 1 < len(le_stations):
        le_r_type, le_r_jd = le_stations[i]
        le_d_type, le_d_jd = le_stations[i + 1]
        se_r_type, se_r_jd = se_stations[i]
        se_d_type, se_d_jd = se_stations[i + 1]

        if le_r_type == "R" and le_d_type == "D":
            le_duration = le_d_jd - le_r_jd
            se_duration = se_d_jd - se_r_jd

            diff_hours = abs(le_duration - se_duration) * 24

            # Retrograde duration should be ~20-24 days
            # LE vs SE should agree within 1 hour
            if 15 < le_duration < 30 and diff_hours < 1.0:
                passed += 1
            elif 15 < le_duration < 30:
                failed += 1
                y, m, d, h = swe.revjul(le_r_jd)
                print(
                    f"  FAIL P3 retro {i // 2} {y}-{m:02d}: LE={le_duration:.2f}d SE={se_duration:.2f}d diff={diff_hours:.2f}h"
                )
            else:
                failed += 1
                print(
                    f"  FAIL P3 retro {i // 2}: duration={le_duration:.1f}d (out of range)"
                )

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Daily Mercury position during retrograde
# ============================================================
print("\n=== P4: Daily position during retrogrades ===")

for i in range(0, min(6, len(le_stations) - 1), 2):  # First 3 retrogrades
    if i + 1 < len(le_stations):
        le_r_type, le_r_jd = le_stations[i]
        le_d_type, le_d_jd = le_stations[i + 1]

        if le_r_type == "R" and le_d_type == "D":
            jd = le_r_jd - 5  # Start 5 days before station R
            end = le_d_jd + 5  # End 5 days after station D
            while jd < end:
                try:
                    le_result = ephem.calc_ut(jd, MERCURY, FLAGS)
                    se_result = swe.calc_ut(jd, MERCURY, FLAGS)

                    le_lon = le_result[0][0]
                    se_lon = se_result[0][0]

                    diff_arcsec = abs(le_lon - se_lon) * 3600
                    if diff_arcsec > 180 * 3600:
                        diff_arcsec = 360 * 3600 - diff_arcsec

                    if diff_arcsec < 0.5:  # Sub-arcsecond during retrograde
                        passed += 1
                    else:
                        failed += 1
                        y, m, d, h = swe.revjul(jd)
                        print(
                            f'  FAIL P4 retro {i // 2} {y}-{m:02d}-{d:02d}: diff={diff_arcsec:.3f}"'
                        )

                except Exception as e:
                    errors += 1
                jd += 1.0

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Mars stations 2020-2030
# ============================================================
print("\n=== P5: Mars station timing 2020-2030 ===")

le_mars = find_stations(MARS, jd_2020, jd_2030, step=1.0)
se_mars = find_stations_se(MARS, jd_2020, jd_2030, step=1.0)

print(f"  LE found {len(le_mars)} Mars stations, SE found {len(se_mars)} stations")

if len(le_mars) == len(se_mars):
    passed += 1
    for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(zip(le_mars, se_mars)):
        if le_type == se_type:
            diff_sec = abs(le_jd - se_jd) * 86400
            if diff_sec < 120.0:  # Mars moves slower, allow 2 min
                passed += 1
            else:
                failed += 1
                y, m, d, h = swe.revjul(le_jd)
                print(
                    f"  FAIL P5 Mars station {i} ({le_type}) {y}-{m:02d}-{d:02d}: diff={diff_sec:.1f}s"
                )
        else:
            failed += 1
else:
    failed += 1
    print(f"  FAIL P5: station count mismatch LE={len(le_mars)} SE={len(se_mars)}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Jupiter stations 2020-2030
# ============================================================
print("\n=== P6: Jupiter station timing 2020-2030 ===")

le_jup = find_stations(JUPITER, jd_2020, jd_2030, step=1.0)
se_jup = find_stations_se(JUPITER, jd_2020, jd_2030, step=1.0)

print(f"  LE found {len(le_jup)} Jupiter stations, SE found {len(se_jup)} stations")

if len(le_jup) == len(se_jup):
    passed += 1
    for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(zip(le_jup, se_jup)):
        if le_type == se_type:
            diff_sec = abs(le_jd - se_jd) * 86400
            le_result = ephem.calc_ut(le_jd, JUPITER, FLAGS)
            se_result = swe.calc_ut(se_jd, JUPITER, FLAGS)
            pos_diff = abs(le_result[0][0] - se_result[0][0]) * 3600

            if diff_sec < 300.0:  # Jupiter very slow near station, allow 5 min
                passed += 1
            else:
                failed += 1
                y, m, d, h = swe.revjul(le_jd)
                print(
                    f'  FAIL P6 Jup station {i} ({le_type}) {y}-{m:02d}-{d:02d}: diff={diff_sec:.1f}s pos_diff={pos_diff:.3f}"'
                )
        else:
            failed += 1
else:
    failed += 1
    print(f"  FAIL P6: station count mismatch LE={len(le_jup)} SE={len(se_jup)}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Saturn stations 2020-2030
# ============================================================
print("\n=== P7: Saturn station timing 2020-2030 ===")

le_sat = find_stations(SATURN, jd_2020, jd_2030, step=1.0)
se_sat = find_stations_se(SATURN, jd_2020, jd_2030, step=1.0)

print(f"  LE found {len(le_sat)} Saturn stations, SE found {len(se_sat)} stations")

if len(le_sat) == len(se_sat):
    passed += 1
    for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(zip(le_sat, se_sat)):
        if le_type == se_type:
            diff_sec = abs(le_jd - se_jd) * 86400
            if diff_sec < 600.0:  # Saturn even slower, allow 10 min
                passed += 1
            else:
                failed += 1
                y, m, d, h = swe.revjul(le_jd)
                print(
                    f"  FAIL P7 Sat station {i} ({le_type}) {y}-{m:02d}-{d:02d}: diff={diff_sec:.1f}s"
                )
        else:
            failed += 1
else:
    failed += 1
    print(f"  FAIL P7: station count mismatch LE={len(le_sat)} SE={len(se_sat)}")

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Mercury speed reversal precision
# ============================================================
print("\n=== P8: Speed at station points (should be ~0) ===")

for i, (le_type, le_jd) in enumerate(le_stations[:10]):  # First 10 stations
    try:
        le_result = ephem.calc_ut(le_jd, MERCURY, FLAGS)
        le_speed = le_result[0][3]

        # Speed at station should be very close to 0
        if abs(le_speed) < 0.001:  # <0.001 deg/day ≈ 3.6"/day
            passed += 1
        else:
            failed += 1
            y, m, d, h = swe.revjul(le_jd)
            print(
                f"  FAIL P8 station {i} ({le_type}) {y}-{m:02d}: speed={le_speed:.6f}°/day"
            )

    except Exception as e:
        errors += 1

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
print(
    f"ROUND 56 FINAL: {passed}/{passed + failed} passed ({100 * passed / (passed + failed):.1f}%)"
)
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
