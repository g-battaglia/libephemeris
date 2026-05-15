#!/usr/bin/env python3
"""
Round 32: XYZ Coordinate Output (FLG_XYZ)
=============================================

Tests Cartesian (X,Y,Z) coordinate output for all major bodies.
FLG_XYZ returns rectangular coordinates instead of spherical (lon,lat,dist).

Parts:
  P1: XYZ ecliptic — all planets, 3 epochs
  P2: XYZ + EQUATORIAL (equatorial rectangular)
  P3: XYZ + J2000
  P4: XYZ + HELCTR (heliocentric XYZ)
  P5: XYZ + BARYCTR (barycentric XYZ)
  P6: XYZ speeds (dx/dy/dz)
  P7: XYZ consistency: verify sqrt(x²+y²+z²) = distance from spherical
"""

from __future__ import annotations

import math
import os
import sys
import time
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris.constants import *

_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)

FLG_BARYCTR = 4

BODIES = [
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
    (CHIRON, "Chiron"),
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanLilith"),
]

HELIO_BODIES = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (CHIRON, "Chiron"),
]

EPOCHS = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 3, 20, 15.5, "2024eq"),
    (1990, 7, 15, 6.0, "1990"),
]

# XYZ tolerance in AU
TOL_XYZ = 1e-5  # AU (~1500 km)
TOL_XYZ_SPEED = 1e-5  # AU/day


class R:
    def __init__(self, name):
        self.name = name
        self.passed = self.failed = self.skipped = 0
        self.failures = []
        self.max_diff = 0.0
        self.max_label = ""

    def ok(self, diff=0.0, label=""):
        self.passed += 1
        if diff > self.max_diff:
            self.max_diff = diff
            self.max_label = label

    def fail(self, msg):
        self.failed += 1
        self.failures.append(msg)

    def skip(self, msg=""):
        self.skipped += 1

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'=' * 70}")
        print(f"  {self.name}: {self.passed}/{total} PASSED ({self.skipped} skip)")
        if self.max_diff > 0:
            print(f"  Max diff: {self.max_diff:.2e} AU ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def compare_xyz(r, se_result, le_result, label, tol=TOL_XYZ):
    se_pos = se_result[0] if isinstance(se_result, tuple) else se_result
    le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

    for i, axis in enumerate(["X", "Y", "Z"]):
        diff = abs(se_pos[i] - le_pos[i])
        if diff > tol:
            r.fail(
                f"{label} {axis}: diff={diff:.2e} AU (SE={se_pos[i]:.10f} LE={le_pos[i]:.10f})"
            )
        else:
            r.ok(diff, f"{label} {axis}")

    # Speed components (indices 3,4,5)
    if len(se_pos) >= 6 and len(le_pos) >= 6:
        for i, axis in enumerate(["dX", "dY", "dZ"]):
            diff = abs(se_pos[i + 3] - le_pos[i + 3])
            if diff > TOL_XYZ_SPEED:
                r.fail(f"{label} {axis}: diff={diff:.2e} AU/d")
            else:
                r.ok(diff, f"{label} {axis}")


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: XYZ ecliptic — all planets, 3 epochs")
    print("=" * 70)

    r = R("P1: XYZ ecliptic")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_XYZ
            label = f"{epoch_name} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_xyz(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: XYZ + EQUATORIAL")
    print("=" * 70)

    r = R("P2: XYZ+EQ")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_XYZ | FLG_EQUATORIAL
            label = f"{epoch_name} {body_name} EQ"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_xyz(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: XYZ + J2000")
    print("=" * 70)

    r = R("P3: XYZ+J2000")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_XYZ | FLG_J2000
            label = f"{epoch_name} {body_name} J2K"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_xyz(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: XYZ + HELCTR")
    print("=" * 70)

    r = R("P4: XYZ+HELCTR")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in HELIO_BODIES:
            flags = FLG_SPEED | FLG_XYZ | FLG_HELCTR
            label = f"{epoch_name} {body_name} HELIO"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_xyz(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: XYZ + BARYCTR")
    print("=" * 70)

    r = R("P5: XYZ+BARY")

    bary_bodies = [
        (b, n) for b, n in BODIES if b not in (MEAN_NODE, TRUE_NODE, MEAN_APOG)
    ]

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in bary_bodies:
            flags = FLG_SPEED | FLG_XYZ | FLG_BARYCTR
            label = f"{epoch_name} {body_name} BARY"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_xyz(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: XYZ speed validation (finite difference)")
    print("=" * 70)

    r = R("P6: XYZ speed FD")

    jd = swe.julday(2024, 3, 20, 15.5)
    dt = 1.0 / 86400.0  # 1 second

    check_bodies = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    for body_id, body_name in check_bodies:
        flags = FLG_SPEED | FLG_XYZ
        label = f"{body_name} XYZ"

        try:
            le_r = ephem.calc_ut(jd, body_id, flags)
            le_before = ephem.calc_ut(jd - dt, body_id, FLG_XYZ)
            le_after = ephem.calc_ut(jd + dt, body_id, FLG_XYZ)

            for i, axis in enumerate(["dX", "dY", "dZ"]):
                reported = le_r[0][i + 3]
                fd = (le_after[0][i] - le_before[0][i]) / (2 * dt)
                fd_diff = abs(reported - fd)
                if fd_diff > 1e-4:
                    r.fail(
                        f"{label} {axis} FD: rep={reported:.10f} fd={fd:.10f} diff={fd_diff:.2e}"
                    )
                else:
                    r.ok(fd_diff, f"{label} {axis} FD")
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: XYZ consistency — sqrt(x²+y²+z²) = dist from spherical")
    print("=" * 70)

    r = R("P7: XYZ consistency")

    jd = swe.julday(2024, 3, 20, 15.5)

    check_bodies = [
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
        (CHIRON, "Chiron"),
    ]

    for body_id, body_name in check_bodies:
        label = f"{body_name}"

        try:
            # Spherical
            sph = ephem.calc_ut(jd, body_id, FLG_SPEED)
            # XYZ
            xyz = ephem.calc_ut(jd, body_id, FLG_SPEED | FLG_XYZ)

            dist_sph = sph[0][2]
            dist_xyz = math.sqrt(xyz[0][0] ** 2 + xyz[0][1] ** 2 + xyz[0][2] ** 2)
            diff = abs(dist_sph - dist_xyz)

            if diff > 1e-10:
                r.fail(
                    f"{label}: sph_dist={dist_sph:.10f} xyz_dist={dist_xyz:.10f} diff={diff:.2e}"
                )
            else:
                r.ok(diff, f"{label}")

            # Also verify lon/lat from XYZ matches spherical
            x, y, z = xyz[0][0], xyz[0][1], xyz[0][2]
            lon_xyz = math.degrees(math.atan2(y, x)) % 360
            lat_xyz = math.degrees(math.atan2(z, math.sqrt(x**2 + y**2)))
            lon_sph = sph[0][0]
            lat_sph = sph[0][1]

            lon_err = abs(lon_xyz - lon_sph)
            if lon_err > 180:
                lon_err = 360 - lon_err
            lat_err = abs(lat_xyz - lat_sph)

            if lon_err > 1e-8 or lat_err > 1e-8:
                r.fail(f"{label} sph↔xyz: Δlon={lon_err:.2e}° Δlat={lat_err:.2e}°")
            else:
                r.ok(lon_err, f"{label} sph↔xyz")

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 32: XYZ Coordinate Output (FLG_XYZ)")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

    for pname, pfn in [
        ("P1", run_part1),
        ("P2", run_part2),
        ("P3", run_part3),
        ("P4", run_part4),
        ("P5", run_part5),
        ("P6", run_part6),
        ("P7", run_part7),
    ]:
        try:
            ok, res = pfn()
            all_results.append((pname, res))
            if not ok:
                all_ok = False
        except Exception as e:
            print(f"\n  {pname} CRASHED: {e}")
            traceback.print_exc()
            all_ok = False

    elapsed = time.time() - start
    print("\n" + "=" * 70)
    print("ROUND 32 FINAL SUMMARY")
    print("=" * 70)
    tp = tf = ts = 0
    for pn, res in all_results:
        st = "PASS" if res.failed == 0 else "FAIL"
        t = res.passed + res.failed
        print(f"  {pn} {res.name}: {res.passed}/{t} ({res.skipped} skip) [{st}]")
        tp += res.passed
        tf += res.failed
        ts += res.skipped
    print(f"\n  TOTAL: {tp}/{tp + tf} PASSED, {tf} FAILED, {ts} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")
    print(f"\n  >>> ROUND 32: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
