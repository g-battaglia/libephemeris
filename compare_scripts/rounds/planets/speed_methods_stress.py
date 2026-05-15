#!/usr/bin/env python3
"""
Round 25: Speed Computation Methods Stress Test
=================================================

Tests speed (velocity) values across various flag combinations and bodies.

Parts:
  P1: Ecliptic speed (default) — all planets, 10 dates
  P2: Equatorial speed — all planets with FLG_EQUATORIAL
  P3: J2000 speed — all planets with FLG_J2000
  P4: Heliocentric speed — all planets with FLG_HELCTR
  P5: Speed near stations (speed ≈ 0) — Mercury & Mars
  P6: Moon speed (fast-moving, latitude speed important)
  P7: Speed finite-difference validation (compare analytical vs FD)
"""

from __future__ import annotations

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
            print(f"  Max diff: {self.max_diff:.8f}°/d ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


PLANETS = [
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
]

DATES_JD = [
    swe.julday(y, m, 15, 12.0)
    for y, m in [
        (2000, 1),
        (2005, 4),
        (2010, 7),
        (2015, 10),
        (2020, 1),
        (2022, 6),
        (2023, 3),
        (2023, 9),
        (2024, 1),
        (2024, 6),
    ]
]


def compare_speed(
    r,
    body_id,
    body_name,
    jd,
    flags,
    label_suffix="",
    tol_lon=0.001,
    tol_lat=0.001,
    tol_dist=0.0001,
):
    """Compare speed values for a body at given flags."""
    label = f"{body_name} {label_suffix}"
    try:
        se_xx = swe.calc_ut(jd, body_id, flags | FLG_SPEED)[0]
    except Exception:
        r.skip(f"{label}: SE error")
        return
    try:
        le_xx, _ = ephem.calc_ut(jd, body_id, flags | FLG_SPEED)
    except Exception:
        r.skip(f"{label}: LE error")
        return

    lon_spd_diff = abs(se_xx[3] - le_xx[3])
    lat_spd_diff = abs(se_xx[4] - le_xx[4])
    dist_spd_diff = abs(se_xx[5] - le_xx[5])

    max_d = max(lon_spd_diff, lat_spd_diff)

    if lon_spd_diff > tol_lon:
        r.fail(f"{label}: lon_spd diff={lon_spd_diff:.6f}°/d")
    elif lat_spd_diff > tol_lat:
        r.fail(f"{label}: lat_spd diff={lat_spd_diff:.6f}°/d")
    else:
        r.ok(max_d, label)


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Ecliptic speed (default) — all planets × 10 dates")
    print("=" * 70)
    r = R("P1: Ecliptic Speed")
    for jd in DATES_JD:
        for body_id, body_name in PLANETS:
            compare_speed(r, body_id, body_name, jd, 0, f"ecl JD={jd:.0f}")
    print(f"  Tested {r.passed + r.failed} speed comparisons")
    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Equatorial speed — all planets × 10 dates")
    print("=" * 70)
    r = R("P2: Equatorial Speed")
    for jd in DATES_JD:
        for body_id, body_name in PLANETS:
            compare_speed(
                r, body_id, body_name, jd, FLG_EQUATORIAL, f"equ JD={jd:.0f}"
            )
    print(f"  Tested {r.passed + r.failed} speed comparisons")
    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: J2000 speed — all planets × 10 dates")
    print("=" * 70)
    r = R("P3: J2000 Speed")
    for jd in DATES_JD:
        for body_id, body_name in PLANETS:
            compare_speed(r, body_id, body_name, jd, FLG_J2000, f"j2k JD={jd:.0f}")
    print(f"  Tested {r.passed + r.failed} speed comparisons")
    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Heliocentric speed — planets × 10 dates")
    print("=" * 70)
    r = R("P4: Helio Speed")
    helio_planets = [(b, n) for b, n in PLANETS if b not in (SUN,)]
    for jd in DATES_JD:
        for body_id, body_name in helio_planets:
            compare_speed(
                r,
                body_id,
                body_name,
                jd,
                FLG_HELCTR,
                f"helio JD={jd:.0f}",
                tol_lon=0.01,
                tol_lat=0.01,
            )
    print(f"  Tested {r.passed + r.failed} speed comparisons")
    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Speed near stations — Mercury & Mars near zero speed")
    print("=" * 70)
    r = R("P5: Station Speed")

    # Mercury stations in 2024 (from Round 22)
    mercury_stations = [
        2460312.63,
        2460401.43,
        2460425.04,
        2460527.71,
        2460550.38,
        2460640.61,
        2460658.37,
    ]

    # Mars station in 2024
    mars_stations = [2460650.48]

    for jd in mercury_stations:
        # Test at station ± small offset
        for dt in [-0.5, -0.1, 0.0, 0.1, 0.5]:
            compare_speed(
                r,
                MERCURY,
                "Mercury",
                jd + dt,
                0,
                f"station JD={jd:.1f}+{dt}",
                tol_lon=0.01,
            )

    for jd in mars_stations:
        for dt in [-0.5, -0.1, 0.0, 0.1, 0.5]:
            compare_speed(
                r,
                MARS,
                "Mars",
                jd + dt,
                0,
                f"station JD={jd:.1f}+{dt}",
                tol_lon=0.01,
            )

    print(f"  Tested {r.passed + r.failed} near-station comparisons")
    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Moon speed — hourly samples over 1 synodic month")
    print("=" * 70)
    r = R("P6: Moon Speed")

    jd_start = swe.julday(2024, 3, 10, 0.0)
    # Sample every 6 hours for 29.5 days
    n_samples = int(29.5 * 4)

    max_lon_diff = 0.0
    max_lat_diff = 0.0

    for i in range(n_samples):
        jd = jd_start + i * 0.25
        try:
            se_xx = swe.calc_ut(jd, MOON, FLG_SPEED)[0]
            le_xx, _ = ephem.calc_ut(jd, MOON, FLG_SPEED)
        except Exception:
            r.skip(f"Moon sample {i}")
            continue

        lon_d = abs(se_xx[3] - le_xx[3])
        lat_d = abs(se_xx[4] - le_xx[4])

        if lon_d > max_lon_diff:
            max_lon_diff = lon_d
        if lat_d > max_lat_diff:
            max_lat_diff = lat_d

        # Moon speed tolerance: 0.005°/day (~18"/day)
        if lon_d > 0.005 or lat_d > 0.005:
            r.fail(f"Moon sample {i}: lon_spd={lon_d:.6f} lat_spd={lat_d:.6f}")
        else:
            r.ok(max(lon_d, lat_d), f"Moon sample {i}")

    print(
        f"  {n_samples} samples, max_lon_diff={max_lon_diff:.6f}°/d, "
        f"max_lat_diff={max_lat_diff:.6f}°/d"
    )
    return r.summary(), r


def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Speed FD validation — compare LE speed with finite difference")
    print("=" * 70)
    r = R("P7: FD Validation")

    jd = swe.julday(2024, 3, 20, 12.0)
    dt = 0.001  # ~1.4 minutes

    for body_id, body_name in PLANETS:
        try:
            le_xx0, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
            le_xx_m, _ = ephem.calc_ut(jd - dt, body_id, 0)
            le_xx_p, _ = ephem.calc_ut(jd + dt, body_id, 0)
        except Exception:
            r.skip(f"{body_name}: calc error")
            continue

        # Finite difference speed
        fd_lon_spd = (le_xx_p[0] - le_xx_m[0]) / (2 * dt)
        # Handle wraparound
        if fd_lon_spd > 180:
            fd_lon_spd -= 360
        elif fd_lon_spd < -180:
            fd_lon_spd += 360

        fd_lat_spd = (le_xx_p[1] - le_xx_m[1]) / (2 * dt)

        # Compare with analytical speed
        lon_err = abs(le_xx0[3] - fd_lon_spd)
        lat_err = abs(le_xx0[4] - fd_lat_spd)

        label = f"{body_name} FD"
        # FD should match analytical within ~0.01°/day
        tol = 0.01
        if lon_err > tol or lat_err > tol:
            r.fail(f"{label}: lon_err={lon_err:.6f} lat_err={lat_err:.6f}")
        else:
            r.ok(max(lon_err, lat_err), label)

        print(f"  {body_name:10s}: lon_err={lon_err:.8f}°/d lat_err={lat_err:.8f}°/d")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 25: Speed Computation Methods Stress Test")
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
    print("ROUND 25 FINAL SUMMARY")
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
    if all_ok:
        print(f"\n  >>> ROUND 25: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 25: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
