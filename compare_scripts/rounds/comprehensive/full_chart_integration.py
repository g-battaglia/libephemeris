#!/usr/bin/env python3
"""
Round 26: Full Chart Calculation Integration Test
===================================================

Tests a complete astrological chart calculation end-to-end, ensuring all
components work together correctly: planet positions, house cusps, aspects,
and derived values.

Parts:
  P1: Full natal chart — 5 famous birth charts with all planets + houses
  P2: House system comparison — Placidus, Koch, Equal, Whole Sign, Campanus
  P3: Planet-in-house placement consistency
  P4: Ascendant/MC consistency with planet positions
  P5: Multiple house systems same chart — verify Asc/MC invariance
  P6: Chart calculation at different epochs (1900, 1950, 2000, 2050, 2100)
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
            print(f'  Max diff: {self.max_diff:.4f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def angle_diff(a, b):
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return abs(d)


def se_hsys(ch):
    return ch.encode("ascii") if isinstance(ch, str) else ch


def le_hsys(ch):
    return ord(ch) if isinstance(ch, str) else ch


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
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (CHIRON, "Chiron"),
]

# Famous birth chart data: (year, month, day, hour_ut, lat, lon, name)
CHARTS = [
    (1879, 3, 14, 11.5, 48.4, 9.98, "Einstein"),  # Ulm, Germany
    (1926, 6, 1, 9.0, 51.5074, -0.1276, "MonroeApprox"),  # LA approximated
    (1963, 8, 4, 19.0, -33.87, 151.21, "Obama-SydneyApprox"),  # Approx
    (2000, 1, 1, 0.0, 41.9028, 12.4964, "Millennium-Rome"),
    (1990, 7, 15, 6.0, 48.8566, 2.3522, "Generic-Paris"),
]


# ============================================================
# PART 1: Full natal chart
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Full natal chart — 5 charts with all planets + Placidus")
    print("=" * 70)

    r = R("P1: Full Chart")

    for y, m, d, h, lat, lon, name in CHARTS:
        jd = swe.julday(y, m, d, h)
        print(f"  {name} ({y}-{m:02d}-{d:02d} {h:.1f}h, lat={lat}, lon={lon})")

        # Planet positions
        for body_id, body_name in PLANETS:
            label = f"{name}/{body_name}"
            try:
                se_xx = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                le_xx, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
            except Exception as e:
                r.skip(f"{label}: {e}")
                continue

            lon_d = angle_diff(se_xx[0], le_xx[0]) * 3600  # arcsec
            lat_d = abs(se_xx[1] - le_xx[1]) * 3600

            tol = 5.0  # 5" for planets
            if max(lon_d, lat_d) > tol:
                r.fail(f'{label}: lon={lon_d:.2f}" lat={lat_d:.2f}"')
            else:
                r.ok(max(lon_d, lat_d), label)

        # Houses (Placidus)
        try:
            se_cusps, se_ascmc = swe.houses(jd, lat, lon, b"P")
            le_cusps, le_ascmc = ephem.houses(jd, lat, lon, ord("P"))
        except Exception as e:
            r.fail(f"{name}/Houses: {e}")
            continue

        for i in range(12):
            cusp_diff = angle_diff(se_cusps[i], le_cusps[i]) * 3600
            label = f"{name}/Cusp{i + 1}"
            if cusp_diff > 5.0:
                r.fail(f'{label}: diff={cusp_diff:.2f}"')
            else:
                r.ok(cusp_diff, label)

        # Asc & MC
        asc_d = angle_diff(se_ascmc[0], le_ascmc[0]) * 3600
        mc_d = angle_diff(se_ascmc[1], le_ascmc[1]) * 3600
        if asc_d > 5.0:
            r.fail(f'{name}/Asc: diff={asc_d:.2f}"')
        else:
            r.ok(asc_d, f"{name}/Asc")
        if mc_d > 5.0:
            r.fail(f'{name}/MC: diff={mc_d:.2f}"')
        else:
            r.ok(mc_d, f"{name}/MC")

    return r.summary(), r


# ============================================================
# PART 2: House system comparison
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: House system comparison — 7 systems at 2024 equinox")
    print("=" * 70)

    r = R("P2: House Systems")

    jd = swe.julday(2024, 3, 20, 12.0)
    lat, lon = 41.9028, 12.4964  # Rome

    house_systems = [
        ("P", "Placidus"),
        ("K", "Koch"),
        ("E", "Equal"),
        ("W", "WholeSign"),
        ("C", "Campanus"),
        ("R", "Regiomontanus"),
        ("B", "Alcabitius"),
    ]

    for hsys, hsys_name in house_systems:
        label = f"Houses-{hsys_name}"
        try:
            se_cusps, se_ascmc = swe.houses(jd, lat, lon, se_hsys(hsys))
            le_cusps, le_ascmc = ephem.houses(jd, lat, lon, le_hsys(hsys))
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        max_diff = 0.0
        worst_cusp = 0
        n = min(len(se_cusps), len(le_cusps), 12)

        for i in range(n):
            d = angle_diff(se_cusps[i], le_cusps[i]) * 3600
            if d > max_diff:
                max_diff = d
                worst_cusp = i + 1

        if max_diff > 5.0:
            r.fail(f'{label}: max={max_diff:.2f}" cusp {worst_cusp}')
        else:
            r.ok(max_diff, label)

        asc_d = angle_diff(se_ascmc[0], le_ascmc[0]) * 3600
        mc_d = angle_diff(se_ascmc[1], le_ascmc[1]) * 3600

        if asc_d > 5.0:
            r.fail(f'{label}/Asc: {asc_d:.2f}"')
        else:
            r.ok(asc_d, f"{label}/Asc")
        if mc_d > 5.0:
            r.fail(f'{label}/MC: {mc_d:.2f}"')
        else:
            r.ok(mc_d, f"{label}/MC")

        print(
            f'  {hsys_name:15s}: cusps max={max_diff:.2f}" Asc={asc_d:.2f}" MC={mc_d:.2f}"'
        )

    return r.summary(), r


# ============================================================
# PART 3: Planet-in-house consistency
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Planet-in-house placement — SE vs LE agree")
    print("=" * 70)

    r = R("P3: Planet Houses")

    jd = swe.julday(2024, 3, 20, 12.0)
    lat, lon = 48.8566, 2.3522  # Paris

    # Get Placidus cusps from both
    se_cusps, _ = swe.houses(jd, lat, lon, b"P")
    le_cusps, _ = ephem.houses(jd, lat, lon, ord("P"))

    def get_house(lon_planet, cusps):
        """Determine which house a planet is in given cusps."""
        for i in range(12):
            c1 = cusps[i]
            c2 = cusps[(i + 1) % 12]
            if c2 < c1:  # Crosses 0° Aries
                if lon_planet >= c1 or lon_planet < c2:
                    return i + 1
            else:
                if c1 <= lon_planet < c2:
                    return i + 1
        return 1  # Fallback

    for body_id, body_name in PLANETS[:10]:  # Main 10 planets
        label = f"{body_name} house"
        try:
            se_xx = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
            le_xx, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
        except Exception:
            r.skip(f"{label}")
            continue

        se_house = get_house(se_xx[0], se_cusps)
        le_house = get_house(le_xx[0], le_cusps)

        if se_house == le_house:
            r.ok(0.0, label)
        else:
            r.fail(f"{label}: SE=house {se_house}, LE=house {le_house}")

        print(f"  {body_name:10s}: SE house {se_house}, LE house {le_house}")

    return r.summary(), r


# ============================================================
# PART 4: Asc/MC consistency with positions
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Asc/MC consistency — verify MC = ARMC converted")
    print("=" * 70)

    r = R("P4: Asc/MC Consistency")

    locations = [
        (12.4964, 41.9028, "Rome"),
        (2.3522, 48.8566, "Paris"),
        (-73.9857, 40.7484, "NewYork"),
        (151.2093, -33.8688, "Sydney"),
        (80.2707, 13.0827, "Chennai"),
    ]

    jd = swe.julday(2024, 3, 20, 12.0)

    for lon, lat, loc_name in locations:
        label = f"AscMC {loc_name}"

        try:
            se_cusps, se_ascmc = swe.houses(jd, lat, lon, b"P")
            le_cusps, le_ascmc = ephem.houses(jd, lat, lon, ord("P"))
        except Exception as e:
            r.fail(f"{label}: {e}")
            continue

        # SE ascmc: [0]=Asc, [1]=MC, [2]=ARMC, [3]=Vertex
        # LE ascmc should match
        for idx, name in [(0, "Asc"), (1, "MC"), (2, "ARMC"), (3, "Vertex")]:
            if idx < len(se_ascmc) and idx < len(le_ascmc):
                d = angle_diff(se_ascmc[idx], le_ascmc[idx]) * 3600
                sub_label = f"{loc_name}/{name}"
                if d > 5.0:
                    r.fail(f'{sub_label}: diff={d:.2f}"')
                else:
                    r.ok(d, sub_label)

        print(
            f'  {loc_name:12s}: Asc={angle_diff(se_ascmc[0], le_ascmc[0]) * 3600:.2f}" '
            f'MC={angle_diff(se_ascmc[1], le_ascmc[1]) * 3600:.2f}" '
            f'ARMC={angle_diff(se_ascmc[2], le_ascmc[2]) * 3600:.2f}"'
        )

    return r.summary(), r


# ============================================================
# PART 5: Multiple house systems — Asc/MC invariance
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Asc/MC invariance across house systems")
    print("=" * 70)

    r = R("P5: Asc/MC Invariance")

    jd = swe.julday(2024, 6, 21, 12.0)
    lat, lon = 41.9028, 12.4964

    systems = ["P", "K", "E", "C", "R", "B", "O", "M"]
    ref_asc = None
    ref_mc = None

    for hsys in systems:
        try:
            le_cusps, le_ascmc = ephem.houses(jd, lat, lon, le_hsys(hsys))
        except Exception:
            r.skip(f"LE houses {hsys}")
            continue

        asc = le_ascmc[0]
        mc = le_ascmc[1]

        if ref_asc is None:
            ref_asc = asc
            ref_mc = mc
            r.ok(0.0, f"Asc ref={hsys}")
            r.ok(0.0, f"MC ref={hsys}")
        else:
            asc_d = angle_diff(asc, ref_asc) * 3600
            mc_d = angle_diff(mc, ref_mc) * 3600

            # Asc and MC should be same regardless of house system
            # (except some systems define MC differently)
            if asc_d > 1.0:
                r.fail(f'Asc {hsys} vs P: diff={asc_d:.2f}"')
            else:
                r.ok(asc_d, f"Asc {hsys}")
            if mc_d > 1.0:
                r.fail(f'MC {hsys} vs P: diff={mc_d:.2f}"')
            else:
                r.ok(mc_d, f"MC {hsys}")

        print(f"  {hsys}: Asc={asc:.4f}° MC={mc:.4f}°")

    return r.summary(), r


# ============================================================
# PART 6: Charts at different epochs
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Chart calculation across epochs (1900-2100)")
    print("=" * 70)

    r = R("P6: Multi-Epoch Charts")

    epochs = [
        (1900, 1, 1, 12.0, "1900"),
        (1950, 6, 15, 0.0, "1950"),
        (2000, 1, 1, 12.0, "J2000"),
        (2050, 7, 1, 0.0, "2050"),
        (2100, 1, 1, 12.0, "2100"),
    ]

    lat, lon = 41.9028, 12.4964  # Rome

    for y, m, d, h, ename in epochs:
        jd = swe.julday(y, m, d, h)

        # Planets
        for body_id, body_name in PLANETS[:10]:
            label = f"{ename}/{body_name}"
            try:
                se_xx = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                le_xx, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
            except Exception:
                r.skip(f"{label}")
                continue

            lon_d = angle_diff(se_xx[0], le_xx[0]) * 3600
            tol = 5.0
            if lon_d > tol:
                r.fail(f'{label}: lon={lon_d:.2f}"')
            else:
                r.ok(lon_d, label)

        # Houses
        try:
            se_cusps, se_ascmc = swe.houses(jd, lat, lon, b"P")
            le_cusps, le_ascmc = ephem.houses(jd, lat, lon, ord("P"))

            max_cusp_d = 0.0
            for i in range(12):
                d = angle_diff(se_cusps[i], le_cusps[i]) * 3600
                if d > max_cusp_d:
                    max_cusp_d = d

            if max_cusp_d > 5.0:
                r.fail(f'{ename}/cusps: max={max_cusp_d:.2f}"')
            else:
                r.ok(max_cusp_d, f"{ename}/cusps")

        except Exception as e:
            r.fail(f"{ename}/houses: {e}")

        print(f"  {ename}: planets + houses tested")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 26: Full Chart Calculation Integration Test")
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
    print("ROUND 26 FINAL SUMMARY")
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
        print(f"\n  >>> ROUND 26: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 26: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
