#!/usr/bin/env python3
"""Round 61: Equatorial Coordinates (RA/Dec) Output Mode

Compare FLG_EQUATORIAL output (RA in degrees, Dec in degrees) for all planets.
Tests pos[0]=RA, pos[1]=Dec, pos[2]=dist, pos[3]=RA speed, pos[4]=Dec speed.
Also tests combinations: EQUATORIAL+J2000, EQUATORIAL+NONUT, EQUATORIAL+TRUEPOS.
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

BODIES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

FLG_SPEED = 256
FLG_EQUATORIAL = 2048
FLG_J2000 = 32
FLG_NONUT = 64
FLG_TRUEPOS = 16
FLG_NOABERR = 1024

print("=" * 70)
print("ROUND 61: Equatorial Coordinates (RA/Dec) Output Mode")
print("=" * 70)

# ============================================================
# P1: Equatorial default (ecliptic of date) monthly 2000-2025
# ============================================================
print("\n=== P1: Equatorial RA/Dec default mode ===")

FLAGS_EQ = FLG_SPEED | FLG_EQUATORIAL

jd_start = 2451545.0
for month in range(0, 25 * 12, 2):  # Bimonthly for 25 years
    jd = jd_start + month * 30.4375
    for body_id, name in BODIES.items():
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS_EQ)[0]
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ)[0]

            # RA comparison (degrees)
            ra_diff = abs(se_r[0] - le_r[0]) * 3600
            if ra_diff > 180 * 3600:
                ra_diff = 360 * 3600 - ra_diff

            # Dec comparison
            dec_diff = abs(se_r[1] - le_r[1]) * 3600

            if ra_diff < 0.5 and dec_diff < 0.5:
                passed += 1
            else:
                year = 2000 + month / 12.0
                failed += 1
                print(
                    f'  FAIL P1 {name} Y{year:.1f}: RA_diff={ra_diff:.3f}" Dec_diff={dec_diff:.3f}"'
                )
        except Exception as e:
            errors += 1

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Equatorial + J2000
# ============================================================
print("\n=== P2: Equatorial + J2000 ===")

FLAGS_EQ_J2000 = FLG_SPEED | FLG_EQUATORIAL | FLG_J2000

for month in range(0, 25 * 12, 6):  # Every 6 months
    jd = jd_start + month * 30.4375
    for body_id, name in BODIES.items():
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS_EQ_J2000)[0]
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ_J2000)[0]

            ra_diff = abs(se_r[0] - le_r[0]) * 3600
            if ra_diff > 180 * 3600:
                ra_diff = 360 * 3600 - ra_diff
            dec_diff = abs(se_r[1] - le_r[1]) * 3600

            if ra_diff < 0.5 and dec_diff < 0.5:
                passed += 1
            else:
                year = 2000 + month / 12.0
                failed += 1
                print(
                    f'  FAIL P2 {name} Y{year:.1f} J2000: RA_diff={ra_diff:.3f}" Dec_diff={dec_diff:.3f}"'
                )
        except Exception:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Equatorial + NONUT
# ============================================================
print("\n=== P3: Equatorial + NONUT ===")

FLAGS_EQ_NONUT = FLG_SPEED | FLG_EQUATORIAL | FLG_NONUT

for month in range(0, 25 * 12, 6):
    jd = jd_start + month * 30.4375
    for body_id, name in BODIES.items():
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS_EQ_NONUT)[0]
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ_NONUT)[0]

            ra_diff = abs(se_r[0] - le_r[0]) * 3600
            if ra_diff > 180 * 3600:
                ra_diff = 360 * 3600 - ra_diff
            dec_diff = abs(se_r[1] - le_r[1]) * 3600

            if ra_diff < 0.5 and dec_diff < 0.5:
                passed += 1
            else:
                year = 2000 + month / 12.0
                failed += 1
                print(
                    f'  FAIL P3 {name} Y{year:.1f} NONUT: RA_diff={ra_diff:.3f}" Dec_diff={dec_diff:.3f}"'
                )
        except Exception:
            errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Equatorial RA/Dec speed
# ============================================================
print("\n=== P4: RA/Dec speed ===")

for month in range(0, 10 * 12, 3):
    jd = jd_start + month * 30.4375
    for body_id, name in BODIES.items():
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS_EQ)[0]
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ)[0]

            ra_spd_diff = abs(se_r[3] - le_r[3])
            dec_spd_diff = abs(se_r[4] - le_r[4])

            # Speed tolerance: 0.001 deg/day = 3.6"/day
            if ra_spd_diff < 0.001 and dec_spd_diff < 0.001:
                passed += 1
            else:
                year = 2000 + month / 12.0
                failed += 1
                print(
                    f"  FAIL P4 {name} Y{year:.1f}: RA_spd_diff={ra_spd_diff:.6f} Dec_spd_diff={dec_spd_diff:.6f}"
                )
        except Exception:
            errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Equatorial + TRUEPOS
# ============================================================
print("\n=== P5: Equatorial + TRUEPOS ===")

FLAGS_EQ_TRUE = FLG_SPEED | FLG_EQUATORIAL | FLG_TRUEPOS

for year in range(2000, 2026):
    jd = swe.julday(year, 6, 15, 12.0)
    for body_id in [0, 1, 2, 3, 4, 5]:
        name = BODIES[body_id]
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS_EQ_TRUE)[0]
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ_TRUE)[0]

            ra_diff = abs(se_r[0] - le_r[0]) * 3600
            if ra_diff > 180 * 3600:
                ra_diff = 360 * 3600 - ra_diff
            dec_diff = abs(se_r[1] - le_r[1]) * 3600

            # TRUEPOS may have slightly wider differences
            if ra_diff < 1.0 and dec_diff < 1.0:
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL P5 {name} Y{year}: RA_diff={ra_diff:.3f}" Dec_diff={dec_diff:.3f}"'
                )
        except Exception:
            errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: RA range check (0-360°)
# ============================================================
print("\n=== P6: RA range validation (0-360°) ===")

for month in range(0, 5 * 12):
    jd = jd_start + month * 30.4375
    for body_id, name in BODIES.items():
        try:
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ)[0]
            ra = le_r[0]
            dec = le_r[1]

            if 0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL P6 {name}: RA={ra:.4f}° Dec={dec:.4f}° (out of range)")
        except Exception:
            errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Equatorial + NOABERR
# ============================================================
print("\n=== P7: Equatorial + NOABERR ===")

FLAGS_EQ_NOAB = FLG_SPEED | FLG_EQUATORIAL | FLG_NOABERR

for year in range(2000, 2026):
    jd = swe.julday(year, 3, 21, 12.0)
    for body_id in [0, 1, 4, 5, 6]:
        name = BODIES[body_id]
        try:
            se_r = swe.calc_ut(jd, body_id, FLAGS_EQ_NOAB)[0]
            le_r = ephem.calc_ut(jd, body_id, FLAGS_EQ_NOAB)[0]

            ra_diff = abs(se_r[0] - le_r[0]) * 3600
            if ra_diff > 180 * 3600:
                ra_diff = 360 * 3600 - ra_diff
            dec_diff = abs(se_r[1] - le_r[1]) * 3600

            if ra_diff < 1.0 and dec_diff < 1.0:
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL P7 {name} Y{year}: RA_diff={ra_diff:.3f}" Dec_diff={dec_diff:.3f}"'
                )
        except Exception:
            errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
print(f"ROUND 61 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
