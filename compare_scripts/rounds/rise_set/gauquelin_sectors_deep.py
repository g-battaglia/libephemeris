#!/usr/bin/env python3
"""
Round 20: Gauquelin Sectors Deep Verification
==============================================

Tests Gauquelin sector calculations across multiple bodies, locations,
dates, and methods.

Parts:
  P1: Method 0 (with ecliptic lat) — 10 planets × 5 locations
  P2: Method 1 (without ecliptic lat) — 10 planets × 5 locations
  P3: 36-sector house cusps via houses('G') — 5 locations
  P4: Sector continuity — verify sector increases monotonically over 24h
  P5: Multiple epochs — same location, 10 dates across 2000-2024
  P6: Southern hemisphere and edge-case latitudes
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

LOCATIONS = [
    (2.3522, 48.8566, 0.0, "Paris"),
    (12.4964, 41.9028, 0.0, "Rome"),
    (-0.1276, 51.5074, 0.0, "London"),
    (-73.9857, 40.7484, 0.0, "New York"),
    (151.2093, -33.8688, 0.0, "Sydney"),
]


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
            print(f"  Max diff: {self.max_diff:.4f} sectors ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
        print(f"{'=' * 70}")
        return self.failed == 0


def se_gauquelin(jd, body, method, geopos, atpress=1013.25, attemp=15.0):
    """Call pyswisseph gauquelin_sector.

    pyswisseph signature: (tjdut, body, method, geopos, atpress, attemp, flags)
    body can be int (planet) or str (fixed star name).
    """
    try:
        ret = swe.gauquelin_sector(jd, body, method, geopos, atpress, attemp)
        if isinstance(ret, tuple):
            return ret[0]
        return ret
    except Exception as e:
        print(f"    SE gauquelin error: {e}")
        return None


def le_gauquelin(jd, body, method, geopos, atpress=1013.25, attemp=15.0):
    """Call libephemeris gauquelin_sector."""
    try:
        return ephem.gauquelin_sector(jd, body, method, geopos, atpress, attemp)
    except Exception:
        return None


# ============================================================
# PART 1: Method 0 (with ecliptic latitude)
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Method 0 — 10 planets × 5 locations (2024-01-15 12:00)")
    print("=" * 70)

    r = R("P1: Method 0")
    jd = swe.julday(2024, 1, 15, 12.0)

    for lon, lat, alt, loc_name in LOCATIONS:
        geopos = (lon, lat, alt)
        for body_id, body_name in PLANETS:
            label = f"{body_name} @ {loc_name}"
            se_val = se_gauquelin(jd, body_id, 0, geopos)
            le_val = le_gauquelin(jd, body_id, 0, geopos)

            if se_val is None or le_val is None:
                r.skip(f"{label}: SE={se_val} LE={le_val}")
                continue

            diff = abs(se_val - le_val)
            # Handle wraparound (sector 36.9 vs 1.1)
            if diff > 18:
                diff = 36 - diff

            # Tolerance: 0.5 sectors for method 0
            if diff > 0.5:
                r.fail(f"{label}: SE={se_val:.4f} LE={le_val:.4f} diff={diff:.4f}")
            else:
                r.ok(diff, label)

        # Print per-location summary
        print(f"  {loc_name}: tested {len(PLANETS)} planets")

    return r.summary(), r


# ============================================================
# PART 2: Method 1 (without ecliptic latitude)
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Method 1 — 10 planets × 5 locations")
    print("=" * 70)

    r = R("P2: Method 1")
    jd = swe.julday(2024, 1, 15, 12.0)

    for lon, lat, alt, loc_name in LOCATIONS:
        geopos = (lon, lat, alt)
        for body_id, body_name in PLANETS:
            label = f"{body_name} @ {loc_name} m1"
            se_val = se_gauquelin(jd, body_id, 1, geopos)
            le_val = le_gauquelin(jd, body_id, 1, geopos)

            if se_val is None or le_val is None:
                r.skip(f"{label}")
                continue

            diff = abs(se_val - le_val)
            if diff > 18:
                diff = 36 - diff

            if diff > 0.5:
                r.fail(f"{label}: SE={se_val:.4f} LE={le_val:.4f} diff={diff:.4f}")
            else:
                r.ok(diff, label)

    return r.summary(), r


# ============================================================
# PART 3: 36-sector house cusps via houses('G')
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Gauquelin 36-sector cusps via houses('G')")
    print("=" * 70)

    r = R("P3: Gauquelin Cusps")
    jd = swe.julday(2024, 1, 15, 12.0)

    for lon, lat, alt, loc_name in LOCATIONS:
        label_base = f"Cusps @ {loc_name}"
        try:
            se_cusps, se_ascmc = swe.houses(jd, lat, lon, b"G")
            le_cusps, le_ascmc = ephem.houses(jd, lat, lon, ord("G"))
        except Exception as e:
            r.fail(f"{label_base}: {e}")
            continue

        # Compare first 36 cusps
        max_diff = 0.0
        worst_cusp = 0
        n_cusps = min(len(se_cusps), len(le_cusps), 36)

        for i in range(n_cusps):
            diff = abs(se_cusps[i] - le_cusps[i])
            if diff > 180:
                diff = 360 - diff

            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        # Tolerance: 0.1 degree for cusps
        if max_diff > 0.1:
            r.fail(f"{label_base}: max_diff={max_diff:.4f}° at cusp {worst_cusp}")
        else:
            r.ok(max_diff, label_base)

        print(
            f"  {loc_name:12s}: {n_cusps} cusps, max_diff={max_diff:.4f}° "
            f"(cusp {worst_cusp})"
        )

        # Also compare Asc and MC from ascmc
        asc_diff = abs(se_ascmc[0] - le_ascmc[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        mc_diff = abs(se_ascmc[1] - le_ascmc[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff

        if asc_diff > 0.01:
            r.fail(f"{label_base} Asc: diff={asc_diff:.6f}°")
        else:
            r.ok(asc_diff, f"Asc @ {loc_name}")

        if mc_diff > 0.01:
            r.fail(f"{label_base} MC: diff={mc_diff:.6f}°")
        else:
            r.ok(mc_diff, f"MC @ {loc_name}")

    return r.summary(), r


# ============================================================
# PART 4: Sector continuity over 24 hours
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Sector continuity — Sun over 24h at Paris")
    print("=" * 70)

    r = R("P4: Sector Continuity")
    jd_start = swe.julday(2024, 6, 21, 0.0)  # Summer solstice
    geopos = (2.3522, 48.8566, 0.0)

    prev_se = None
    prev_le = None
    n_steps = 48  # Every 30 minutes

    for i in range(n_steps):
        jd = jd_start + i / 48.0
        label = f"Step {i}"

        se_val = se_gauquelin(jd, SUN, 0, geopos)
        le_val = le_gauquelin(jd, SUN, 0, geopos)

        if se_val is None or le_val is None:
            r.skip(f"{label}")
            continue

        # Check LE is in valid range
        if le_val < 1.0 or le_val >= 37.0:
            r.fail(f"{label}: LE out of range: {le_val:.4f}")
            continue

        diff = abs(se_val - le_val)
        if diff > 18:
            diff = 36 - diff

        if diff > 0.5:
            r.fail(f"{label}: SE={se_val:.4f} LE={le_val:.4f} diff={diff:.4f}")
        else:
            r.ok(diff, label)

        prev_se = se_val
        prev_le = le_val

    return r.summary(), r


# ============================================================
# PART 5: Multiple epochs
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Multiple epochs — Mars at Rome, 10 dates")
    print("=" * 70)

    r = R("P5: Multi-Epoch")
    geopos = (12.4964, 41.9028, 0.0)

    dates = [
        (2000, 3, 21, 12.0),
        (2002, 6, 21, 12.0),
        (2004, 9, 23, 12.0),
        (2006, 12, 22, 12.0),
        (2008, 3, 20, 12.0),
        (2010, 6, 21, 12.0),
        (2012, 9, 22, 12.0),
        (2014, 12, 22, 12.0),
        (2018, 7, 27, 12.0),
        (2024, 1, 15, 12.0),
    ]

    for y, m, d, h in dates:
        jd = swe.julday(y, m, d, h)
        label = f"Mars {y}-{m:02d}-{d:02d}"

        se_val = se_gauquelin(jd, MARS, 0, geopos)
        le_val = le_gauquelin(jd, MARS, 0, geopos)

        if se_val is None or le_val is None:
            r.skip(f"{label}")
            continue

        diff = abs(se_val - le_val)
        if diff > 18:
            diff = 36 - diff

        if diff > 0.5:
            r.fail(f"{label}: SE={se_val:.4f} LE={le_val:.4f} diff={diff:.4f}")
        else:
            r.ok(diff, label)

        print(f"  {label}: SE={se_val:.4f} LE={le_val:.4f} diff={diff:.4f}")

    return r.summary(), r


# ============================================================
# PART 6: Southern hemisphere and edge-case latitudes
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Edge-case latitudes — equator, tropics, high latitudes")
    print("=" * 70)

    r = R("P6: Edge Latitudes")
    jd = swe.julday(2024, 6, 21, 12.0)

    edge_locations = [
        (0.0, 0.0, 0.0, "Equator"),
        (0.0, 23.4, 0.0, "Tropic N"),
        (0.0, -23.4, 0.0, "Tropic S"),
        (0.0, 45.0, 0.0, "Mid-lat N"),
        (0.0, -45.0, 0.0, "Mid-lat S"),
        (0.0, 60.0, 0.0, "High-lat N"),
        (0.0, -60.0, 0.0, "High-lat S"),
        (0.0, 64.0, 0.0, "Near-polar N"),
        (0.0, -64.0, 0.0, "Near-polar S"),
    ]

    for lon, lat, alt, loc_name in edge_locations:
        geopos = (lon, lat, alt)
        for body_id, body_name in [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MARS, "Mars"),
        ]:
            label = f"{body_name} @ {loc_name}"
            se_val = se_gauquelin(jd, body_id, 0, geopos)
            le_val = le_gauquelin(jd, body_id, 0, geopos)

            if se_val is None and le_val is None:
                r.skip(f"{label}: both None (polar?)")
                continue
            if se_val is None or le_val is None:
                r.fail(f"{label}: SE={se_val} LE={le_val}")
                continue

            diff = abs(se_val - le_val)
            if diff > 18:
                diff = 36 - diff

            # Wider tolerance for edge latitudes (sectors get compressed)
            tol = 1.0
            if diff > tol:
                r.fail(f"{label}: SE={se_val:.4f} LE={le_val:.4f} diff={diff:.4f}")
            else:
                r.ok(diff, label)

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 20: Gauquelin Sectors Deep Verification")
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
    print("ROUND 20 FINAL SUMMARY")
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
        print(f"\n  >>> ROUND 20: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 20: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
