#!/usr/bin/env python3
"""Round 85: XYZ Output Consistency

Tests FLG_XYZ flag which outputs cartesian coordinates (X, Y, Z) instead of
(lon, lat, dist). Verifies all planets across multiple dates and flag combos.
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

FLG_SWIEPH = 2
FLG_SPEED = 256
FLG_HELCTR = 8
FLG_J2000 = 32
FLG_NONUT = 64
FLG_EQUATORIAL = 2048
FLG_XYZ = 4096

print("=" * 70)
print("ROUND 85: XYZ Output Consistency")
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


def check_xyz(label, se, le, tol_au=0.0001):
    global passed, failed
    for i, axis in enumerate(["X", "Y", "Z"]):
        diff = abs(se[0][i] - le[0][i])
        if diff < tol_au:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL {label} {axis}: SE={se[0][i]:.8f} LE={le[0][i]:.8f} diff={diff:.8f} AU"
            )


# ============================================================
# P1: Geocentric ecliptic XYZ
# ============================================================
print("\n=== P1: Geocentric ecliptic XYZ ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_XYZ
for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            check_xyz(label, se, le)
        except Exception as e:
            errors += 1
            if errors <= 3:
                print(f"  ERROR {label}: {e}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Heliocentric XYZ
# ============================================================
print("\n=== P2: Heliocentric XYZ ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_XYZ | FLG_HELCTR
for body_id, name in bodies[1:]:  # Skip Sun
    for year, jd in test_jds:
        label = f"{name} helio {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            check_xyz(label, se, le)
        except Exception as e:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: J2000 XYZ
# ============================================================
print("\n=== P3: J2000 ecliptic XYZ ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_XYZ | FLG_J2000 | FLG_NONUT
for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} J2000 {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            check_xyz(label, se, le)
        except Exception as e:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Equatorial XYZ
# ============================================================
print("\n=== P4: Equatorial XYZ ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_XYZ | FLG_EQUATORIAL
for body_id, name in [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
]:
    for year, jd in test_jds:
        label = f"{name} eq {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            check_xyz(label, se, le)
        except Exception as e:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: XYZ speeds (velocity components)
# ============================================================
print("\n=== P5: XYZ velocity components ===")
flags = FLG_SWIEPH | FLG_SPEED | FLG_XYZ
for body_id, name in [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
]:
    for year, jd in test_jds:
        label = f"{name} vel {year}"
        try:
            se = swe.calc_ut(jd, body_id, flags)
            le = ephem.calc_ut(jd, body_id, flags)
            for i, axis in enumerate(["vX", "vY", "vZ"]):
                diff = abs(se[0][3 + i] - le[0][3 + i])
                if diff < 0.0001:  # AU/day
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL {label} {axis}: SE={se[0][3 + i]:.8f} LE={le[0][3 + i]:.8f}"
                    )
        except Exception as e:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: XYZ self-consistency — verify X²+Y²+Z² = dist²
# ============================================================
print("\n=== P6: XYZ distance self-consistency ===")
for body_id, name in [(swe.SUN, "Sun"), (swe.MOON, "Moon"), (swe.MARS, "Mars")]:
    for year, jd in test_jds[:10]:
        label = f"{name} dist_check {year}"
        try:
            # Get spherical
            le_sph = ephem.calc_ut(jd, body_id, FLG_SWIEPH | FLG_SPEED)
            dist_sph = le_sph[0][2]
            # Get XYZ
            le_xyz = ephem.calc_ut(
                jd, body_id, FLG_SWIEPH | FLG_SPEED | FLG_XYZ
            )
            dist_xyz = math.sqrt(
                le_xyz[0][0] ** 2 + le_xyz[0][1] ** 2 + le_xyz[0][2] ** 2
            )
            ratio = dist_xyz / dist_sph if dist_sph != 0 else 999
            if abs(ratio - 1.0) < 0.0001:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL {label}: sph_dist={dist_sph:.8f} xyz_dist={dist_xyz:.8f} ratio={ratio:.8f}"
                )
        except Exception as e:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 85 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
