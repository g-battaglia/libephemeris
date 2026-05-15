#!/usr/bin/env python3
"""
Round 33: Radians Output Mode (FLG_RADIANS)
===============================================

Tests radians output mode for all major bodies across multiple flag combos.
FLG_RADIANS converts lon/lat/speed output from degrees to radians.

Parts:
  P1: Radians ecliptic — all planets, 3 epochs
  P2: Radians + EQUATORIAL
  P3: Radians + J2000
  P4: Radians + HELCTR
  P5: Radians consistency: verify rad = deg * pi/180
  P6: Radians speeds consistency
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
]

EPOCHS = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 3, 20, 15.5, "2024eq"),
    (1990, 7, 15, 6.0, "1990"),
]

# Radians tolerance: 2" in radians = 2/(3600*180/pi) ≈ 9.7e-6 rad
TOL_RAD = 2.0 / 3600.0 * math.pi / 180.0  # ~9.7e-6 rad
TOL_DIST = 1e-5  # AU
TOL_SPEED_RAD = 0.001 * math.pi / 180.0  # 0.001 deg/day in rad/day
TOL_MEANLILITH_LAT_RAD = 25.0 / 3600.0 * math.pi / 180.0


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
            arcsec = self.max_diff * 180 / math.pi * 3600
            print(
                f'  Max diff: {self.max_diff:.2e} rad ({arcsec:.3f}") ({self.max_label})'
            )
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def compare_rad(r, se_result, le_result, label, body_id):
    se_pos = se_result[0] if isinstance(se_result, tuple) else se_result
    le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

    # Longitude (radians)
    lon_diff = abs(se_pos[0] - le_pos[0])
    if lon_diff > math.pi:
        lon_diff = 2 * math.pi - lon_diff

    # Latitude
    lat_diff = abs(se_pos[1] - le_pos[1])
    lat_tol = TOL_MEANLILITH_LAT_RAD if body_id == MEAN_APOG else TOL_RAD

    # Distance
    dist_diff = abs(se_pos[2] - le_pos[2])

    if lon_diff > TOL_RAD:
        arcsec = lon_diff * 180 / math.pi * 3600
        r.fail(f'{label} lon: {arcsec:.3f}" ({lon_diff:.2e} rad)')
    elif lat_diff > lat_tol:
        arcsec = lat_diff * 180 / math.pi * 3600
        r.fail(f'{label} lat: {arcsec:.3f}" ({lat_diff:.2e} rad)')
    elif dist_diff > TOL_DIST:
        r.fail(f"{label} dist: {dist_diff:.8f} AU")
    else:
        r.ok(lon_diff, f"{label}")

    # Speeds
    if len(se_pos) >= 4 and len(le_pos) >= 4:
        spd_diff = abs(se_pos[3] - le_pos[3])
        if spd_diff > TOL_SPEED_RAD:
            r.fail(f"{label} spd: {spd_diff:.2e} rad/d")
        else:
            r.ok(spd_diff, f"{label} spd")


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Radians ecliptic — all planets, 3 epochs")
    print("=" * 70)

    r = R("P1: Radians ecliptic")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_RADIANS
            label = f"{epoch_name} {body_name} RAD"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_rad(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Radians + EQUATORIAL")
    print("=" * 70)

    r = R("P2: Radians+EQ")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_RADIANS | FLG_EQUATORIAL
            label = f"{epoch_name} {body_name} RAD+EQ"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_rad(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Radians + J2000")
    print("=" * 70)

    r = R("P3: Radians+J2000")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_RADIANS | FLG_J2000
            label = f"{epoch_name} {body_name} RAD+J2K"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_rad(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Radians + HELCTR")
    print("=" * 70)

    r = R("P4: Radians+HELCTR")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in HELIO_BODIES:
            flags = FLG_SPEED | FLG_RADIANS | FLG_HELCTR
            label = f"{epoch_name} {body_name} RAD+HELIO"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_rad(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Radians consistency — verify rad = deg × π/180")
    print("=" * 70)

    r = R("P5: Radians consistency")

    jd = swe.julday(2024, 3, 20, 15.5)

    for body_id, body_name in BODIES:
        label = f"{body_name}"

        try:
            deg_r = ephem.calc_ut(jd, body_id, FLG_SPEED)
            rad_r = ephem.calc_ut(jd, body_id, FLG_SPEED | FLG_RADIANS)

            # lon: deg → rad
            expected_lon_rad = deg_r[0][0] * math.pi / 180.0
            actual_lon_rad = rad_r[0][0]
            diff = abs(expected_lon_rad - actual_lon_rad)
            if diff > 1e-12:
                r.fail(
                    f"{label} lon: exp={expected_lon_rad:.12f} got={actual_lon_rad:.12f} diff={diff:.2e}"
                )
            else:
                r.ok(diff, f"{label} lon")

            # lat: deg → rad
            expected_lat_rad = deg_r[0][1] * math.pi / 180.0
            actual_lat_rad = rad_r[0][1]
            diff = abs(expected_lat_rad - actual_lat_rad)
            if diff > 1e-12:
                r.fail(f"{label} lat: diff={diff:.2e}")
            else:
                r.ok(diff, f"{label} lat")

            # dist: should be identical
            diff = abs(deg_r[0][2] - rad_r[0][2])
            if diff > 1e-15:
                r.fail(f"{label} dist: diff={diff:.2e}")
            else:
                r.ok(diff, f"{label} dist")

            # speed: deg/day → rad/day
            expected_spd_rad = deg_r[0][3] * math.pi / 180.0
            actual_spd_rad = rad_r[0][3]
            diff = abs(expected_spd_rad - actual_spd_rad)
            if diff > 1e-12:
                r.fail(f"{label} spd: diff={diff:.2e}")
            else:
                r.ok(diff, f"{label} spd")

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Radians speed consistency — all 6 components")
    print("=" * 70)

    r = R("P6: Radians speed full")

    jd = swe.julday(2024, 3, 20, 15.5)

    for body_id, body_name in BODIES:
        label = f"{body_name}"

        try:
            deg_r = ephem.calc_ut(jd, body_id, FLG_SPEED)
            rad_r = ephem.calc_ut(jd, body_id, FLG_SPEED | FLG_RADIANS)

            # Speed components 3,4 should be deg→rad, 5 (dist speed) unchanged
            for i, comp in enumerate(["lon_spd", "lat_spd", "dist_spd"]):
                deg_val = deg_r[0][i + 3]
                rad_val = rad_r[0][i + 3]

                if i < 2:  # angular speed
                    expected = deg_val * math.pi / 180.0
                else:  # distance speed (AU/day)
                    expected = deg_val

                diff = abs(expected - rad_val)
                if diff > 1e-12:
                    r.fail(
                        f"{label} {comp}: exp={expected:.12f} got={rad_val:.12f} diff={diff:.2e}"
                    )
                else:
                    r.ok(diff, f"{label} {comp}")

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 33: Radians Output Mode (FLG_RADIANS)")
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
    print("ROUND 33 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 33: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
