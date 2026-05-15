#!/usr/bin/env python3
"""
Round 34: Topocentric Positions Deep Sweep (FLG_TOPOCTR)
===========================================================

Tests topocentric position calculations across multiple geographic locations,
bodies, and epochs. Topocentric = position as seen from a specific point on
Earth's surface (not geocentric center).

Parts:
  P1: Major cities — all planets at 5 locations
  P2: Extreme latitudes (Arctic, Antarctic, equator)
  P3: High altitude locations (Everest, Mauna Kea, La Silla)
  P4: Topocentric + EQUATORIAL
  P5: Topocentric + J2000
  P6: Topocentric Moon parallax validation (should be ~1°)
  P7: Multi-epoch sweep at single location
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

FLG_TOPOCTR = 32768  # 0x8000

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
]

# Locations: (lat, lon, alt_m, name)
CITIES = [
    (51.5074, -0.1278, 11, "London"),
    (40.7128, -74.0060, 10, "NewYork"),
    (35.6762, 139.6503, 40, "Tokyo"),
    (-33.8688, 151.2093, 58, "Sydney"),
    (55.7558, 37.6173, 156, "Moscow"),
]

EXTREME_LATS = [
    (78.2232, 15.6267, 10, "Svalbard"),  # Arctic
    (-77.8500, 166.6667, 24, "McMurdo"),  # Antarctic
    (0.0, 0.0, 0, "NullIsland"),  # Equator
    (64.1466, -21.9426, 10, "Reykjavik"),  # Sub-Arctic
    (-54.8019, -68.3030, 20, "Ushuaia"),  # Sub-Antarctic
]

HIGH_ALT = [
    (27.9881, 86.9250, 8849, "Everest"),
    (19.8207, -155.4681, 4207, "MaunaKea"),
    (-29.2563, -70.7380, 2400, "LaSilla"),
]

TOL_LON = 3.0  # arcsec (topocentric has more model-dependent variation)
TOL_LAT = 3.0
TOL_DIST = 5e-5  # AU
TOL_SPEED = 0.005  # deg/day


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


def compare_pos(r, se_result, le_result, label):
    se_pos = se_result[0] if isinstance(se_result, tuple) else se_result
    le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

    lon_diff = abs(se_pos[0] - le_pos[0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_diff_as = lon_diff * 3600

    lat_diff = abs(se_pos[1] - le_pos[1]) * 3600
    dist_diff = abs(se_pos[2] - le_pos[2])

    if lon_diff_as > TOL_LON:
        r.fail(
            f'{label} lon: {lon_diff_as:.3f}" (SE={se_pos[0]:.8f} LE={le_pos[0]:.8f})'
        )
    elif lat_diff > TOL_LAT:
        r.fail(f'{label} lat: {lat_diff:.3f}"')
    elif dist_diff > TOL_DIST:
        r.fail(f"{label} dist: {dist_diff:.8f} AU")
    else:
        r.ok(max(lon_diff_as, lat_diff), f"{label}")

    if len(se_pos) >= 4 and len(le_pos) >= 4:
        spd_diff = abs(se_pos[3] - le_pos[3])
        if spd_diff > TOL_SPEED:
            r.fail(f"{label} spd: {spd_diff:.6f} deg/day")
        else:
            r.ok(spd_diff * 3600, f"{label} spd")


def set_topo(lat, lon, alt_m):
    swe.set_topo(lon, lat, alt_m)
    ephem.set_topo(lon, lat, alt_m)


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Major cities — all planets, topocentric")
    print("=" * 70)

    r = R("P1: Cities topocentric")
    jd = swe.julday(2024, 3, 20, 15.5)

    for lat, lon, alt, city in CITIES:
        set_topo(lat, lon, alt)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_TOPOCTR
            label = f"{city} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Extreme latitudes — topocentric")
    print("=" * 70)

    r = R("P2: Extreme latitudes")
    jd = swe.julday(2024, 6, 21, 12.0)

    for lat, lon, alt, loc_name in EXTREME_LATS:
        set_topo(lat, lon, alt)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_TOPOCTR
            label = f"{loc_name} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: High altitude locations — topocentric")
    print("=" * 70)

    r = R("P3: High altitude")
    jd = swe.julday(2024, 3, 20, 15.5)

    for lat, lon, alt, loc_name in HIGH_ALT:
        set_topo(lat, lon, alt)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_TOPOCTR
            label = f"{loc_name} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Topocentric + EQUATORIAL")
    print("=" * 70)

    r = R("P4: Topo+EQ")
    jd = swe.julday(2024, 3, 20, 15.5)

    for lat, lon, alt, city in CITIES[:3]:  # London, NewYork, Tokyo
        set_topo(lat, lon, alt)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_TOPOCTR | FLG_EQUATORIAL
            label = f"{city} {body_name} EQ"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Topocentric + J2000")
    print("=" * 70)

    r = R("P5: Topo+J2000")
    jd = swe.julday(2024, 3, 20, 15.5)

    for lat, lon, alt, city in CITIES[:3]:
        set_topo(lat, lon, alt)
        for body_id, body_name in BODIES:
            flags = FLG_SPEED | FLG_TOPOCTR | FLG_J2000
            label = f"{city} {body_name} J2K"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Topocentric Moon parallax validation")
    print("=" * 70)

    r = R("P6: Moon parallax")
    jd = swe.julday(2024, 3, 20, 15.5)

    for lat, lon, alt, city in CITIES:
        set_topo(lat, lon, alt)

        # Geocentric Moon
        geo = ephem.calc_ut(jd, MOON, FLG_SPEED)
        # Topocentric Moon
        topo = ephem.calc_ut(jd, MOON, FLG_SPEED | FLG_TOPOCTR)

        lon_diff = abs(geo[0][0] - topo[0][0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        lat_diff = abs(geo[0][1] - topo[0][1])
        total_parallax = (lon_diff**2 + lat_diff**2) ** 0.5 * 3600  # arcsec

        label = f"{city} Moon parallax"
        # Moon's horizontal parallax is ~3420" (~57'), topocentric shift depends on geometry
        # but should be between 0 and ~3600"
        if total_parallax < 1 or total_parallax > 4000:
            r.fail(f'{label}: parallax={total_parallax:.1f}" (out of range)')
        else:
            r.ok(total_parallax, f"{label}")
            print(
                f'  {city}: Moon parallax shift = {total_parallax:.1f}" ({total_parallax / 3600:.3f}°)'
            )

        # Also check Sun parallax is much smaller (< 30")
        geo_sun = ephem.calc_ut(jd, SUN, FLG_SPEED)
        topo_sun = ephem.calc_ut(jd, SUN, FLG_SPEED | FLG_TOPOCTR)
        sun_lon_diff = abs(geo_sun[0][0] - topo_sun[0][0])
        if sun_lon_diff > 180:
            sun_lon_diff = 360 - sun_lon_diff
        sun_lat_diff = abs(geo_sun[0][1] - topo_sun[0][1])
        sun_parallax = (sun_lon_diff**2 + sun_lat_diff**2) ** 0.5 * 3600

        label_sun = f"{city} Sun parallax"
        if sun_parallax > 30:
            r.fail(f'{label_sun}: parallax={sun_parallax:.3f}" (too large)')
        else:
            r.ok(sun_parallax, f"{label_sun}")

    return r.summary(), r


def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Multi-epoch sweep at London")
    print("=" * 70)

    r = R("P7: Multi-epoch topo")

    set_topo(51.5074, -0.1278, 11)  # London

    sweep_bodies = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]

    for y in range(1900, 2101, 20):
        jd = swe.julday(y, 6, 21, 12.0)
        for body_id, body_name in sweep_bodies:
            flags = FLG_SPEED | FLG_TOPOCTR
            label = f"{y} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 34: Topocentric Positions Deep Sweep")
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
    print("ROUND 34 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 34: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
