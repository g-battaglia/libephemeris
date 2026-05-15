#!/usr/bin/env python3
"""Round 86: Radians Output Mode

Tests FLG_RADIANS flag which outputs positions in radians instead of degrees.
Verifies consistency with degree output and across flag combinations.
"""

from __future__ import annotations

import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLG_SWIEPH = 2
FLG_SPEED = 256
FLG_J2000 = 32
FLG_NONUT = 64
FLG_EQUATORIAL = 2048
FLG_RADIANS = 8192

print("=" * 70)
print("ROUND 86: Radians Output Mode")
print("=" * 70)

bodies = [
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
]

test_jds = [(y, swe.julday(y, 1, 15, 12.0)) for y in range(1980, 2030, 2)]

# ============================================================
# P1: Radians SE vs LE comparison
# ============================================================
print("\n=== P1: Radians ecliptic SE vs LE ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_RADIANS
for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            for i, coord in enumerate(["lon", "lat", "dist"]):
                if i < 2:
                    diff = abs(se[0][i] - le[0][i])
                    if diff > math.pi:
                        diff = 2 * math.pi - diff
                    diff_arcsec = math.degrees(diff) * 3600.0
                    if diff_arcsec < 1.0:
                        passed += 1
                    else:
                        failed += 1
                        print(
                            f'  FAIL {label} {coord}: SE={se[0][i]:.8f} LE={le[0][i]:.8f} diff={diff_arcsec:.2f}"'
                        )
                else:
                    if le[0][i] != 0:
                        ratio = se[0][i] / le[0][i]
                        if abs(ratio - 1.0) < 0.0001:
                            passed += 1
                        else:
                            failed += 1
                            print(
                                f"  FAIL {label} dist: SE={se[0][i]:.8f} LE={le[0][i]:.8f}"
                            )
                    else:
                        passed += 1
        except Exception as e:
            errors += 1
            if errors <= 3:
                print(f"  ERROR {label}: {e}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Radians self-consistency (radians = degrees * pi/180)
# ============================================================
print("\n=== P2: Radians = degrees * pi/180 self-consistency ===")
for body_id, name in bodies:
    for year, jd in test_jds[:10]:
        label = f"{name} {year}"
        try:
            le_deg = ephem.calc_ut(jd, body_id, FLG_SWIEPH | FLG_SPEED)
            le_rad = ephem.calc_ut(
                jd, body_id, FLG_SWIEPH | FLG_SPEED | FLG_RADIANS
            )
            for i, coord in enumerate(["lon", "lat"]):
                expected_rad = math.radians(le_deg[0][i])
                diff = abs(le_rad[0][i] - expected_rad)
                if diff > math.pi:
                    diff = 2 * math.pi - diff
                diff_arcsec = math.degrees(diff) * 3600.0
                if diff_arcsec < 0.001:  # sub-milliarcsecond
                    passed += 1
                else:
                    failed += 1
                    print(
                        f'  FAIL {label} {coord}: rad={le_rad[0][i]:.10f} expected={expected_rad:.10f} diff={diff_arcsec:.6f}"'
                    )
            # Distance should be identical
            if abs(le_rad[0][2] - le_deg[0][2]) < 1e-12:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL {label} dist: rad={le_rad[0][2]:.10f} deg={le_deg[0][2]:.10f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Radians speed self-consistency
# ============================================================
print("\n=== P3: Speed in radians ===")
for body_id, name in [(swe.SUN, "Sun"), (swe.MOON, "Moon"), (swe.MARS, "Mars")]:
    for year, jd in test_jds[:10]:
        label = f"{name} {year}"
        try:
            le_deg = ephem.calc_ut(jd, body_id, FLG_SWIEPH | FLG_SPEED)
            le_rad = ephem.calc_ut(
                jd, body_id, FLG_SWIEPH | FLG_SPEED | FLG_RADIANS
            )
            for i, coord in enumerate(["lon_spd", "lat_spd"]):
                expected_rad_spd = math.radians(le_deg[0][3 + i])
                diff = abs(le_rad[0][3 + i] - expected_rad_spd)
                diff_arcsec = math.degrees(diff) * 3600.0
                if diff_arcsec < 0.001:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL {label} {coord}: rad={le_rad[0][3 + i]:.10f} expected={expected_rad_spd:.10f}"
                    )
        except Exception as e:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Radians + equatorial
# ============================================================
print("\n=== P4: Radians + equatorial ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_RADIANS | FLG_EQUATORIAL
for body_id, name in [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
]:
    for year, jd in test_jds:
        label = f"{name} eq {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            for i, coord in enumerate(["RA", "Dec"]):
                diff = abs(se[0][i] - le[0][i])
                if diff > math.pi:
                    diff = 2 * math.pi - diff
                diff_arcsec = math.degrees(diff) * 3600.0
                if diff_arcsec < 1.5:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f'  FAIL {label} {coord}: SE={se[0][i]:.8f} LE={le[0][i]:.8f} diff={diff_arcsec:.2f}"'
                    )
        except Exception as e:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Radians + J2000
# ============================================================
print("\n=== P5: Radians + J2000 ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_RADIANS | FLG_J2000 | FLG_NONUT
for body_id, name in [(swe.SUN, "Sun"), (swe.MOON, "Moon"), (swe.MARS, "Mars")]:
    for year, jd in test_jds:
        label = f"{name} J2000 {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            for i, coord in enumerate(["lon", "lat"]):
                diff = abs(se[0][i] - le[0][i])
                if diff > math.pi:
                    diff = 2 * math.pi - diff
                diff_arcsec = math.degrees(diff) * 3600.0
                if diff_arcsec < 1.0:
                    passed += 1
                else:
                    failed += 1
                    print(f'  FAIL {label} {coord}: diff={diff_arcsec:.2f}"')
        except Exception as e:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: Range validation (lon in [0, 2π), lat in [-π/2, π/2])
# ============================================================
print("\n=== P6: Radians range validation ===")
for body_id, name in bodies:
    for year, jd in test_jds[:5]:
        label = f"{name} {year}"
        try:
            le = ephem.calc_ut(
                jd, body_id, FLG_SWIEPH | FLG_SPEED | FLG_RADIANS
            )
            lon_rad = le[0][0]
            lat_rad = le[0][1]
            if 0 <= lon_rad < 2 * math.pi:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL {label} lon range: {lon_rad:.6f} (expected [0, 2π))")
            if -math.pi / 2 <= lat_rad <= math.pi / 2:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL {label} lat range: {lat_rad:.6f} (expected [-π/2, π/2])")
        except Exception as e:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 86 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
