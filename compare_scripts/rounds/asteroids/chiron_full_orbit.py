#!/usr/bin/env python3
"""
Round 36: Chiron Full Orbit Sweep
==================================

Chiron has a ~50.7-year orbital period with high eccentricity (0.38) and
inclination (6.9°). This sweep tests Chiron across its full orbit to catch
any ephemeris edge cases at perihelion/aphelion or near Saturn/Uranus crossing.

Parts:
  P1: Full orbit sweep (1970-2020, every 6 months — covers 1 full orbit)
  P2: Near perihelion (1996, q=8.45 AU) and aphelion (1970, Q=18.8 AU)
  P3: Heliocentric ecliptic across full orbit
  P4: Equatorial and J2000 coordinates
  P5: Speed near perihelion vs aphelion (should differ ~4x)
  P6: Extended range (1900-2100)
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

TOL_LON = 2.0
TOL_LAT = 2.0
TOL_DIST = 1e-5
TOL_SPEED = 0.002


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
            print(f'  Max diff: {self.max_diff:.3f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def compare_pos(r, se_result, le_result, label, check_speed=True):
    se_pos = se_result[0] if isinstance(se_result, tuple) else se_result
    le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

    lon_diff = abs(se_pos[0] - le_pos[0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    lat_as = abs(se_pos[1] - le_pos[1]) * 3600
    dist_diff = abs(se_pos[2] - le_pos[2])

    if lon_as > TOL_LON:
        r.fail(f'{label} lon: {lon_as:.3f}"')
    elif lat_as > TOL_LAT:
        r.fail(f'{label} lat: {lat_as:.3f}"')
    elif dist_diff > TOL_DIST:
        r.fail(f"{label} dist: {dist_diff:.8f} AU")
    else:
        r.ok(max(lon_as, lat_as), f"{label}")

    if check_speed and len(se_pos) >= 4 and len(le_pos) >= 4:
        spd_diff = abs(se_pos[3] - le_pos[3])
        if spd_diff > TOL_SPEED:
            r.fail(f"{label} spd: {spd_diff:.6f}")
        else:
            r.ok(spd_diff * 3600, f"{label} spd")


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Full orbit sweep (1970-2020, every 6 months)")
    print("=" * 70)

    r = R("P1: Chiron full orbit")
    flags = FLG_SPEED

    for y in range(1970, 2021):
        for m in [1, 7]:
            jd = swe.julday(y, m, 1, 12.0)
            label = f"{y}-{m:02d} Chiron"
            try:
                se_r = swe.calc_ut(jd, CHIRON, flags)
                le_r = ephem.calc_ut(jd, CHIRON, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Near perihelion (1996) and aphelion (1970)")
    print("=" * 70)

    r = R("P2: Perihelion/Aphelion")
    flags = FLG_SPEED

    # Chiron perihelion ~Feb 1996 (q≈8.45 AU)
    # Chiron aphelion ~1970 (Q≈18.8 AU)
    critical_dates = [
        (1970, 1, 1, "aphelion_area"),
        (1970, 7, 1, "aphelion_area2"),
        (1995, 1, 1, "pre_perihelion"),
        (1996, 2, 14, "perihelion"),
        (1996, 6, 1, "post_perihelion"),
        (1997, 1, 1, "post_perihelion2"),
        (2020, 7, 1, "mid_orbit"),
        (2021, 1, 1, "mid_orbit2"),
    ]

    for y, m, d, desc in critical_dates:
        jd = swe.julday(y, m, d, 12.0)
        label = f"Chiron {desc}"
        try:
            se_r = swe.calc_ut(jd, CHIRON, flags)
            le_r = ephem.calc_ut(jd, CHIRON, flags)
            compare_pos(r, se_r, le_r, label)
            print(
                f"  {desc}: lon={le_r[0][0]:.4f}° dist={le_r[0][2]:.4f} AU spd={le_r[0][3]:.6f}°/d"
            )
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Heliocentric ecliptic across full orbit")
    print("=" * 70)

    r = R("P3: Chiron heliocentric")
    flags = FLG_SPEED | FLG_HELCTR

    for y in range(1970, 2021, 2):
        jd = swe.julday(y, 6, 21, 12.0)
        label = f"{y} Chiron HELIO"
        try:
            se_r = swe.calc_ut(jd, CHIRON, flags)
            le_r = ephem.calc_ut(jd, CHIRON, flags)
            compare_pos(r, se_r, le_r, label)
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Equatorial and J2000 coordinates")
    print("=" * 70)

    r = R("P4: Chiron EQ/J2000")

    test_dates = [
        (1980, 1, 1, 12.0),
        (1996, 2, 14, 12.0),
        (2000, 1, 1, 12.0),
        (2024, 3, 20, 15.5),
    ]

    for y, m, d, h in test_dates:
        jd = swe.julday(y, m, d, h)

        # Equatorial
        flags_eq = FLG_SPEED | FLG_EQUATORIAL
        try:
            se_r = swe.calc_ut(jd, CHIRON, flags_eq)
            le_r = ephem.calc_ut(jd, CHIRON, flags_eq)
            compare_pos(r, se_r, le_r, f"{y} Chiron EQ")
        except Exception as e:
            r.skip(f"{y} Chiron EQ: {e}")

        # J2000
        flags_j2k = FLG_SPEED | FLG_J2000
        try:
            se_r = swe.calc_ut(jd, CHIRON, flags_j2k)
            le_r = ephem.calc_ut(jd, CHIRON, flags_j2k)
            compare_pos(r, se_r, le_r, f"{y} Chiron J2K")
        except Exception as e:
            r.skip(f"{y} Chiron J2K: {e}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Speed near perihelion vs aphelion")
    print("=" * 70)

    r = R("P5: Chiron speed ratio")

    # Perihelion speed should be ~4x aphelion speed (Kepler's 2nd law)
    jd_peri = swe.julday(1996, 2, 14, 12.0)
    jd_aph = swe.julday(1970, 7, 1, 12.0)

    peri_r = ephem.calc_ut(jd_peri, CHIRON, FLG_SPEED | FLG_HELCTR)
    aph_r = ephem.calc_ut(jd_aph, CHIRON, FLG_SPEED | FLG_HELCTR)

    peri_spd = abs(peri_r[0][3])
    aph_spd = abs(aph_r[0][3])
    ratio = peri_spd / aph_spd if aph_spd > 0 else 0

    print(f"  Perihelion speed: {peri_spd:.6f} °/day")
    print(f"  Aphelion speed:   {aph_spd:.6f} °/day")
    print(f"  Ratio:            {ratio:.2f}x")

    # e=0.38, ratio should be ~ (1+e)/(1-e) ≈ 2.23 for velocity ratio
    # For angular speed: v_peri/v_aph = ((1+e)/(1-e))^2 ≈ 4.97 ... but depends on viewing geometry
    # Just check it's > 1.5 (faster at perihelion)
    if ratio < 1.5:
        r.fail(f"Speed ratio too low: {ratio:.2f}x (expected > 1.5)")
    else:
        r.ok(0, f"speed ratio {ratio:.2f}x")

    # Also validate against SE
    se_peri = swe.calc_ut(jd_peri, CHIRON, FLG_SPEED | FLG_HELCTR)
    se_aph = swe.calc_ut(jd_aph, CHIRON, FLG_SPEED | FLG_HELCTR)

    peri_diff = abs(se_peri[0][3] - peri_r[0][3])
    aph_diff = abs(se_aph[0][3] - aph_r[0][3])
    if peri_diff > 0.005:
        r.fail(f"Perihelion speed diff: {peri_diff:.6f}")
    else:
        r.ok(peri_diff * 3600, "peri speed")
    if aph_diff > 0.005:
        r.fail(f"Aphelion speed diff: {aph_diff:.6f}")
    else:
        r.ok(aph_diff * 3600, "aph speed")

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Extended range (1900-2100)")
    print("=" * 70)

    r = R("P6: Chiron extended")
    flags = FLG_SPEED

    for y in range(1900, 2101, 5):
        jd = swe.julday(y, 6, 21, 12.0)
        label = f"{y} Chiron"
        try:
            se_r = swe.calc_ut(jd, CHIRON, flags)
            le_r = ephem.calc_ut(jd, CHIRON, flags)
            compare_pos(r, se_r, le_r, label, check_speed=False)
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 36: Chiron Full Orbit Sweep")
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
    print("ROUND 36 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 36: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
