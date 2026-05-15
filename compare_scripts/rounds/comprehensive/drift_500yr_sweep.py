#!/usr/bin/env python3
"""
Round 18: Multi-Epoch 500-Year Drift Sweep
==========================================

Tests all major celestial bodies across 1600-2100 at regular intervals,
checking for systematic positional drift between libephemeris and the
reference implementation.

Parts:
  P1: Main planets (Sun-Pluto) longitude/latitude/distance at 50-year intervals
  P2: Lunar nodes (mean & true) across centuries
  P3: Mean & True Lilith (Mean Apogee, Osculating Apogee) across centuries
  P4: Planet speeds across epochs (FLG_SPEED)
  P5: Equatorial coordinates (FLG_EQUATORIAL) across epochs
  P6: Heliocentric positions (FLG_HELCTR) across epochs
  P7: J2000 ecliptic (FLG_J2000) across epochs — isolate precession from ephemeris
  P8: Dense sweep 1900-2100 at 1-year intervals for Moon (most sensitive body)
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

# Body definitions
MAIN_PLANETS = [
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

NODE_BODIES = [
    (MEAN_NODE, "Mean Node"),
    (TRUE_NODE, "True Node"),
]

LILITH_BODIES = [
    (MEAN_APOG, "Mean Lilith"),
    (OSCU_APOG, "Oscu Lilith"),
]

# Epoch definitions
EPOCHS_50Y = list(range(1600, 2101, 50))  # 1600, 1650, ..., 2100
EPOCHS_DENSE = list(range(1900, 2101, 1))  # 1900, 1901, ..., 2100


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
            print(f'  Max diff: {self.max_diff:.6f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
        print(f"{'=' * 70}")
        return self.failed == 0


def arcsec(deg):
    """Convert degrees to arcseconds."""
    return abs(deg) * 3600.0


def compare_body_at_epoch(
    body_id,
    body_name,
    year,
    flags,
    r,
    tol_lon_arcsec=1.0,
    tol_lat_arcsec=1.0,
    tol_dist_ppm=100.0,
    label_prefix="",
    skip_dist=False,
):
    """Compare a body's position at Jan 1 of the given year."""
    jd = swe.julday(year, 1, 1, 12.0)
    label = f"{label_prefix}{body_name} {year}"

    try:
        se_pos, se_flag = swe.calc_ut(jd, body_id, flags)
        le_pos, le_flag = ephem.calc_ut(jd, body_id, flags)
    except Exception as e:
        r.fail(f"{label}: EXCEPTION {e}")
        return None

    # Longitude difference (handle 360° wrap)
    dlon = se_pos[0] - le_pos[0]
    if dlon > 180:
        dlon -= 360
    elif dlon < -180:
        dlon += 360
    dlon_as = arcsec(dlon)

    # Latitude difference
    dlat = se_pos[1] - le_pos[1]
    dlat_as = arcsec(dlat)

    # Distance difference (parts per million)
    if se_pos[2] != 0:
        ddist_ppm = abs(se_pos[2] - le_pos[2]) / abs(se_pos[2]) * 1e6
    else:
        ddist_ppm = 0.0

    fails = []
    if dlon_as > tol_lon_arcsec:
        fails.append(f'lon={dlon_as:.3f}"')
    if dlat_as > tol_lat_arcsec:
        fails.append(f'lat={dlat_as:.3f}"')
    if not skip_dist and ddist_ppm > tol_dist_ppm:
        fails.append(f"dist={ddist_ppm:.1f}ppm")

    max_as = max(dlon_as, dlat_as)

    if fails:
        r.fail(f"{label}: {'; '.join(fails)}")
    else:
        r.ok(max_as, label)

    return dlon_as, dlat_as, ddist_ppm


# ============================================================
# PART 1: Main planets across 500 years
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Main planets — ecliptic of date, 1600-2100 at 50-yr intervals")
    print("=" * 70)

    r = R("P1: Multi-Epoch Planets")
    flags = FLG_SWIEPH | FLG_SPEED

    # Tolerances vary by body and epoch distance from J2000
    for body_id, body_name in MAIN_PLANETS:
        drifts = []
        for year in EPOCHS_50Y:
            # Allow looser tolerance for distant epochs
            dist_from_j2000 = abs(year - 2000)
            # Base tolerance: 0.5" near J2000, scales up for distant epochs
            tol_lon = 0.5 + dist_from_j2000 * 0.005  # +0.005"/yr
            tol_lat = 0.5 + dist_from_j2000 * 0.005
            tol_dist = 50 + dist_from_j2000 * 0.5  # ppm

            # Moon uses DE440 while SE uses its own lunar theory — drift
            # is ~0.75"/century from J2000, reaching ~3" at 400 years.
            if body_id == MOON:
                tol_lon = 0.5 + dist_from_j2000 * 0.03  # ~3"/century
                tol_lat = 0.5 + dist_from_j2000 * 0.01

            result = compare_body_at_epoch(
                body_id, body_name, year, flags, r, tol_lon, tol_lat, tol_dist
            )
            if result:
                drifts.append((year, result))

        # Print drift summary for this body
        if drifts:
            max_lon = max(d[1][0] for d in drifts)
            max_lat = max(d[1][1] for d in drifts)
            worst_year_lon = max(drifts, key=lambda d: d[1][0])
            print(
                f'  {body_name:10s}: max_lon={max_lon:8.3f}"  '
                f'max_lat={max_lat:8.3f}"  '
                f"worst_year={worst_year_lon[0]}"
            )

    return r.summary(), r


# ============================================================
# PART 2: Lunar nodes across centuries
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Lunar nodes (Mean & True) — 1600-2100")
    print("=" * 70)

    r = R("P2: Lunar Nodes Multi-Epoch")
    flags = FLG_SWIEPH | FLG_SPEED

    for body_id, body_name in NODE_BODIES:
        drifts = []
        for year in EPOCHS_50Y:
            dist_from_j2000 = abs(year - 2000)
            tol_lon = 1.0 + dist_from_j2000 * 0.01
            tol_lat = 1.0 + dist_from_j2000 * 0.01

            result = compare_body_at_epoch(
                body_id,
                body_name,
                year,
                flags,
                r,
                tol_lon,
                tol_lat,
                500.0,
                skip_dist=True,
            )
            if result:
                drifts.append((year, result))

        if drifts:
            max_lon = max(d[1][0] for d in drifts)
            print(f'  {body_name:12s}: max_lon={max_lon:8.3f}"')

    return r.summary(), r


# ============================================================
# PART 3: Mean & Osculating Lilith across centuries
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Lilith (Mean & Osculating Apogee) — 1600-2100")
    print("=" * 70)

    r = R("P3: Lilith Multi-Epoch")
    flags = FLG_SWIEPH | FLG_SPEED

    for body_id, body_name in LILITH_BODIES:
        drifts = []
        for year in EPOCHS_50Y:
            dist_from_j2000 = abs(year - 2000)
            # Lilith has larger inherent differences
            tol_lon = 5.0 + dist_from_j2000 * 0.05
            # Mean Lilith latitude differs ~15-20" due to different
            # computation models (lunar theory vs geometric mean).
            # Osculating Lilith latitude is more consistent.
            tol_lat = 25.0 if body_id == MEAN_APOG else 5.0

            result = compare_body_at_epoch(
                body_id,
                body_name,
                year,
                flags,
                r,
                tol_lon,
                tol_lat,
                1000.0,
                skip_dist=True,
            )
            if result:
                drifts.append((year, result))

        if drifts:
            max_lon = max(d[1][0] for d in drifts)
            print(f'  {body_name:14s}: max_lon={max_lon:8.3f}"')

    return r.summary(), r


# ============================================================
# PART 4: Planet speeds across epochs
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Speed computation across epochs — 1600-2100")
    print("=" * 70)

    r = R("P4: Speed Multi-Epoch")
    flags = FLG_SWIEPH | FLG_SPEED

    for body_id, body_name in MAIN_PLANETS:
        max_speed_diff = 0.0
        worst_year = 0

        for year in EPOCHS_50Y:
            jd = swe.julday(year, 1, 1, 12.0)
            label = f"Speed {body_name} {year}"

            try:
                se_pos, _ = swe.calc_ut(jd, body_id, flags)
                le_pos, _ = ephem.calc_ut(jd, body_id, flags)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            # Speed is in pos[3] (lon speed, deg/day)
            speed_diff = abs(se_pos[3] - le_pos[3])
            speed_diff_as = speed_diff * 3600  # arcsec/day

            dist_from_j2000 = abs(year - 2000)
            tol_speed = 0.5 + dist_from_j2000 * 0.005  # arcsec/day

            if speed_diff_as > tol_speed:
                r.fail(f'{label}: speed_diff={speed_diff_as:.3f}"/day')
            else:
                r.ok(speed_diff_as, label)

            if speed_diff_as > max_speed_diff:
                max_speed_diff = speed_diff_as
                worst_year = year

        print(
            f'  {body_name:10s}: max_speed_diff={max_speed_diff:8.4f}"/day  '
            f"worst_year={worst_year}"
        )

    return r.summary(), r


# ============================================================
# PART 5: Equatorial coordinates across epochs
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Equatorial coords (RA/Dec) — 1700-2100 at 100-yr intervals")
    print("=" * 70)

    r = R("P5: Equatorial Multi-Epoch")
    flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL
    epochs = list(range(1700, 2101, 100))

    for body_id, body_name in MAIN_PLANETS:
        max_ra_diff = 0.0
        for year in epochs:
            jd = swe.julday(year, 1, 1, 12.0)
            label = f"EQ {body_name} {year}"

            try:
                se_pos, _ = swe.calc_ut(jd, body_id, flags)
                le_pos, _ = ephem.calc_ut(jd, body_id, flags)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            # RA difference (handle 360° wrap)
            dra = se_pos[0] - le_pos[0]
            if dra > 180:
                dra -= 360
            elif dra < -180:
                dra += 360
            dra_as = arcsec(dra)

            # Dec difference
            ddec_as = arcsec(se_pos[1] - le_pos[1])

            dist_from_j2000 = abs(year - 2000)
            tol = 1.0 + dist_from_j2000 * 0.01

            fails = []
            if dra_as > tol:
                fails.append(f'RA={dra_as:.3f}"')
            if ddec_as > tol:
                fails.append(f'Dec={ddec_as:.3f}"')

            if fails:
                r.fail(f"{label}: {'; '.join(fails)}")
            else:
                r.ok(max(dra_as, ddec_as), label)

            if dra_as > max_ra_diff:
                max_ra_diff = dra_as

        print(f'  {body_name:10s}: max_RA_diff={max_ra_diff:8.3f}"')

    return r.summary(), r


# ============================================================
# PART 6: Heliocentric positions across epochs
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Heliocentric ecliptic — 1700-2100 at 100-yr intervals")
    print("=" * 70)

    r = R("P6: Heliocentric Multi-Epoch")
    flags = FLG_SWIEPH | FLG_SPEED | FLG_HELCTR
    epochs = list(range(1700, 2101, 100))

    # Skip Sun (heliocentric Sun is always 0) and Moon (not valid heliocentric)
    helio_planets = [p for p in MAIN_PLANETS if p[0] not in (SUN, MOON)]

    for body_id, body_name in helio_planets:
        max_lon_diff = 0.0
        for year in epochs:
            jd = swe.julday(year, 1, 1, 12.0)
            label = f"Helio {body_name} {year}"

            try:
                se_pos, _ = swe.calc_ut(jd, body_id, flags)
                le_pos, _ = ephem.calc_ut(jd, body_id, flags)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            dlon = se_pos[0] - le_pos[0]
            if dlon > 180:
                dlon -= 360
            elif dlon < -180:
                dlon += 360
            dlon_as = arcsec(dlon)
            dlat_as = arcsec(se_pos[1] - le_pos[1])

            dist_from_j2000 = abs(year - 2000)
            # Heliocentric removes geocentric parallax but ephemeris
            # model differences persist. Mercury and Pluto need more margin.
            tol = 0.3 + dist_from_j2000 * 0.005

            fails = []
            if dlon_as > tol:
                fails.append(f'lon={dlon_as:.3f}"')
            if dlat_as > tol:
                fails.append(f'lat={dlat_as:.3f}"')

            if fails:
                r.fail(f"{label}: {'; '.join(fails)}")
            else:
                r.ok(max(dlon_as, dlat_as), label)

            if dlon_as > max_lon_diff:
                max_lon_diff = dlon_as

        print(f'  {body_name:10s}: max_lon_diff={max_lon_diff:8.3f}"')

    return r.summary(), r


# ============================================================
# PART 7: J2000 ecliptic — isolate precession from ephemeris
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: J2000 ecliptic (no precession) — 1700-2100 at 100-yr intervals")
    print("=" * 70)

    r = R("P7: J2000 Multi-Epoch")
    flags = FLG_SWIEPH | FLG_SPEED | FLG_J2000 | FLG_NONUT
    epochs = list(range(1700, 2101, 100))

    for body_id, body_name in MAIN_PLANETS:
        max_lon_diff = 0.0
        for year in epochs:
            jd = swe.julday(year, 1, 1, 12.0)
            label = f"J2000 {body_name} {year}"

            try:
                se_pos, _ = swe.calc_ut(jd, body_id, flags)
                le_pos, _ = ephem.calc_ut(jd, body_id, flags)
            except Exception as e:
                r.fail(f"{label}: {e}")
                continue

            dlon = se_pos[0] - le_pos[0]
            if dlon > 180:
                dlon -= 360
            elif dlon < -180:
                dlon += 360
            dlon_as = arcsec(dlon)
            dlat_as = arcsec(se_pos[1] - le_pos[1])

            dist_from_j2000 = abs(year - 2000)
            # J2000 frame removes precession but ephemeris model
            # differences remain. Moon and Pluto need wider margin.
            if body_id == MOON:
                tol = 0.5 + dist_from_j2000 * 0.02
            elif body_id == PLUTO:
                tol = 0.5 + dist_from_j2000 * 0.005
            else:
                tol = 0.3 + dist_from_j2000 * 0.003

            fails = []
            if dlon_as > tol:
                fails.append(f'lon={dlon_as:.3f}"')
            if dlat_as > tol:
                fails.append(f'lat={dlat_as:.3f}"')

            if fails:
                r.fail(f"{label}: {'; '.join(fails)}")
            else:
                r.ok(max(dlon_as, dlat_as), label)

            if dlon_as > max_lon_diff:
                max_lon_diff = dlon_as

        print(f'  {body_name:10s}: max_lon_diff={max_lon_diff:8.3f}"')

    return r.summary(), r


# ============================================================
# PART 8: Dense Moon sweep 1900-2100 at 1-year intervals
# ============================================================
def run_part8():
    print("\n" + "=" * 70)
    print("PART 8: Moon dense sweep — 1900-2100, yearly")
    print("  Checking for systematic drift patterns")
    print("=" * 70)

    r = R("P8: Moon Dense Sweep")
    flags = FLG_SWIEPH | FLG_SPEED

    lon_diffs = []
    for year in EPOCHS_DENSE:
        jd = swe.julday(year, 1, 1, 12.0)
        label = f"Moon {year}"

        try:
            se_pos, _ = swe.calc_ut(jd, MOON, flags)
            le_pos, _ = ephem.calc_ut(jd, MOON, flags)
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        dlon = se_pos[0] - le_pos[0]
        if dlon > 180:
            dlon -= 360
        elif dlon < -180:
            dlon += 360
        dlon_as = arcsec(dlon)
        dlat_as = arcsec(se_pos[1] - le_pos[1])

        # Moon tolerance scales with distance from J2000 due to
        # DE440 vs SE lunar theory divergence (~0.75"/century linear
        # trend plus ~0.5" oscillation from the 18.6-year nodal cycle).
        dist_from_j2000 = abs(year - 2000)
        tol = 0.5 + dist_from_j2000 * 0.03  # ~3"/century

        fails = []
        if dlon_as > tol:
            fails.append(f'lon={dlon_as:.3f}"')
        if dlat_as > tol:
            fails.append(f'lat={dlat_as:.3f}"')

        if fails:
            r.fail(f"{label}: {'; '.join(fails)}")
        else:
            r.ok(dlon_as, label)

        lon_diffs.append((year, dlon_as, dlat_as))

    # Print statistics
    if lon_diffs:
        max_lon = max(d[1] for d in lon_diffs)
        avg_lon = sum(d[1] for d in lon_diffs) / len(lon_diffs)
        max_lat = max(d[2] for d in lon_diffs)
        avg_lat = sum(d[2] for d in lon_diffs) / len(lon_diffs)
        worst = max(lon_diffs, key=lambda d: d[1])

        print(
            f'  Moon longitude: avg={avg_lon:.4f}"  max={max_lon:.4f}"  worst_year={worst[0]}'
        )
        print(f'  Moon latitude:  avg={avg_lat:.4f}"  max={max_lat:.4f}"')

        # Check for systematic drift (linear regression on lon diffs)
        n = len(lon_diffs)
        x_mean = sum(d[0] for d in lon_diffs) / n
        y_mean = sum(d[1] for d in lon_diffs) / n
        num = sum((d[0] - x_mean) * (d[1] - y_mean) for d in lon_diffs)
        den = sum((d[0] - x_mean) ** 2 for d in lon_diffs)
        if den > 0:
            slope = num / den
            print(f'  Drift rate: {slope * 100:.6f}"/century (linear fit)')
        else:
            print(f"  Drift rate: N/A")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 18: Multi-Epoch 500-Year Drift Sweep")
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
        ("P8", run_part8),
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
    print("ROUND 18 FINAL SUMMARY")
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
        print(f"\n  >>> ROUND 18: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 18: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
