#!/usr/bin/env python3
"""
Round 21: Rise/Set at Extreme Latitudes Deep Verification
==========================================================

Tests rise, transit, and set calculations at challenging latitudes where
circumpolar, never-rise, and grazing conditions occur.

Parts:
  P1: Sun rise/set across latitudes 0°-70° (normal cases)
  P2: Sun rise/set near polar circle (64°-68°) — summer/winter solstice
  P3: Moon rise/set across latitudes (Moon can be ±5° ecliptic lat)
  P4: Planet rise/set at extreme latitudes
  P5: Transit times (should always work) across all latitudes
  P6: Twilight (civil/nautical/astronomical) at extreme latitudes
  P7: True horizon (elevated observer) at various latitudes
"""

from __future__ import annotations

import os
import sys
import time
import traceback
import math

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
            print(f"  Max diff: {self.max_diff:.1f}s ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def jd_diff_seconds(jd1, jd2):
    """Difference in seconds between two JDs."""
    return abs(jd1 - jd2) * 86400.0


def se_rise_trans(jd, body, geopos, rsmi, atpress=1013.25, attemp=15.0):
    """Call pyswisseph rise_trans.

    pyswisseph sig: (tjdut, body, rsmi, geopos, atpress, attemp, flags)
    Returns: (res_flag, tret_tuple) where tret_tuple[0] = JD of event
    res_flag: 0 = found, -2 = circumpolar
    """
    try:
        res, tret = swe.rise_trans(jd, body, rsmi, geopos, atpress, attemp)
        if res == -2:
            return None  # circumpolar
        return tret[0]
    except Exception:
        return None


def le_rise_trans(jd, body, geopos, rsmi, atpress=1013.25, attemp=15.0):
    """Call libephemeris rise_trans.

    libephemeris sig: (tjdut, body, rsmi, geopos, atpress, attemp, flags)
    Returns: (retflag, tret) where retflag=-2 means circumpolar, tret[0] = JD
    """
    try:
        retflag, tret = ephem.rise_trans(
            jd, body, rsmi, list(geopos), atpress, attemp
        )
        if retflag == -2 or tret[0] == 0.0:
            return None  # circumpolar
        return tret[0]
    except Exception:
        return None


def se_rise_trans_true_hor(
    jd, body, geopos, rsmi, horiz_alt=0.0, atpress=1013.25, attemp=15.0
):
    """Call pyswisseph rise_trans_true_hor.

    pyswisseph sig: (tjdut, body, rsmi, geopos, atpress, attemp, horhgt, flags)
    """
    try:
        res, tret = swe.rise_trans_true_hor(
            jd, body, rsmi, geopos, atpress, attemp, horiz_alt
        )
        if res == -2:
            return None
        return tret[0]
    except Exception:
        return None


def le_rise_trans_true_hor(
    jd, body, geopos, rsmi, horiz_alt=0.0, atpress=1013.25, attemp=15.0
):
    """Call libephemeris rise_trans_true_hor.

    libephemeris sig: (tjdut, body, rsmi, geopos, atpress, attemp, horhgt, flags)
    Returns: (retflag, tret) where tret[0] = JD
    """
    try:
        retflag, tret = ephem.rise_trans_true_hor(
            jd, body, rsmi, list(geopos), atpress, attemp, horiz_alt
        )
        if retflag == -2 or tret[0] == 0.0:
            return None
        return tret[0]
    except Exception:
        return None


# RSMI flags (from SE docs)
CALC_RISE = 1
CALC_SET = 2
CALC_MTRANSIT = 4
CALC_ITRANSIT = 8
BIT_DISC_CENTER = 256
BIT_DISC_BOTTOM = 8192
BIT_NO_REFRACTION = 512
BIT_CIVIL_TWILIGHT = 1024
BIT_NAUTIC_TWILIGHT = 2048
BIT_ASTRO_TWILIGHT = 4096


# ============================================================
# PART 1: Sun rise/set across latitudes 0°-70°
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Sun rise/set across latitudes 0°-70° (equinox & solstice)")
    print("=" * 70)

    r = R("P1: Sun Rise/Set Normal")

    # Test at spring equinox and summer solstice
    dates = [
        (2024, 3, 20, 12.0, "Equinox"),
        (2024, 6, 21, 0.0, "Solstice"),
        (2024, 12, 21, 0.0, "Winter Sol"),
    ]

    latitudes = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0, 60.0, 63.0, 65.0, 67.0, 70.0]

    for y, m, d, h, dname in dates:
        jd = swe.julday(y, m, d, h)
        for lat in latitudes:
            geopos = (0.0, lat, 0.0)

            for rsmi, evname in [
                (CALC_RISE, "Rise"),
                (CALC_SET, "Set"),
            ]:
                label = f"Sun {evname} lat={lat}° {dname}"
                se_val = se_rise_trans(jd, SUN, geopos, rsmi)
                le_val = le_rise_trans(jd, SUN, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: both None (circumpolar?)")
                    continue
                if se_val is None or le_val is None:
                    # One returns result, other doesn't — might be edge case
                    r.fail(f"{label}: SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)

                # Tolerance: 30s for normal latitudes, 120s for >60°
                tol = 120.0 if lat > 60.0 else 30.0
                if diff_s > tol:
                    r.fail(
                        f"{label}: diff={diff_s:.1f}s (SE={se_val:.6f} LE={le_val:.6f})"
                    )
                else:
                    r.ok(diff_s, label)

        print(f"  {dname}: tested {len(latitudes)} latitudes")

    return r.summary(), r


# ============================================================
# PART 2: Sun rise/set near polar circle (64°-68°)
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Sun near polar circle — fine latitude sweep")
    print("=" * 70)

    r = R("P2: Sun Polar Circle")

    # Fine sweep near polar circle during summer/winter
    jd_summer = swe.julday(2024, 6, 21, 0.0)
    jd_winter = swe.julday(2024, 12, 21, 0.0)

    # Latitudes near Arctic circle (~66.56°)
    latitudes = [64.0, 64.5, 65.0, 65.5, 66.0, 66.3, 66.5, 66.6, 67.0, 67.5, 68.0]

    for jd, season in [(jd_summer, "Summer"), (jd_winter, "Winter")]:
        for lat in latitudes:
            geopos = (25.0, lat, 0.0)  # Use longitude 25° (Finland area)

            for rsmi, evname in [
                (CALC_RISE, "Rise"),
                (CALC_SET, "Set"),
            ]:
                label = f"Sun {evname} lat={lat}° {season}"
                se_val = se_rise_trans(jd, SUN, geopos, rsmi)
                le_val = le_rise_trans(jd, SUN, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: both circumpolar/never-rise")
                    continue
                if se_val is None or le_val is None:
                    # Disagreement on circumpolar status — important to note
                    r.fail(f"{label}: DISAGREE SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)
                # Wider tolerance near polar circle — grazing conditions
                tol = 300.0
                if diff_s > tol:
                    r.fail(f"{label}: diff={diff_s:.1f}s")
                else:
                    r.ok(diff_s, label)

        print(f"  {season}: tested {len(latitudes)} latitudes")

    # Also test southern polar circle
    jd_june = swe.julday(2024, 6, 21, 0.0)
    for lat in [-64.0, -66.0, -66.5, -67.0, -68.0]:
        geopos = (0.0, lat, 0.0)
        for rsmi, evname in [(CALC_RISE, "Rise"), (CALC_SET, "Set")]:
            label = f"Sun {evname} lat={lat}° S-Winter"
            se_val = se_rise_trans(jd_june, SUN, geopos, rsmi)
            le_val = le_rise_trans(jd_june, SUN, geopos, rsmi)

            if se_val is None and le_val is None:
                r.skip(f"{label}: both circumpolar")
                continue
            if se_val is None or le_val is None:
                r.fail(f"{label}: DISAGREE SE={se_val} LE={le_val}")
                continue

            diff_s = jd_diff_seconds(se_val, le_val)
            tol = 300.0
            if diff_s > tol:
                r.fail(f"{label}: diff={diff_s:.1f}s")
            else:
                r.ok(diff_s, label)

    print(f"  Southern: tested 5 latitudes")

    return r.summary(), r


# ============================================================
# PART 3: Moon rise/set across latitudes
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Moon rise/set across latitudes (5 dates)")
    print("=" * 70)

    r = R("P3: Moon Rise/Set")

    dates = [
        (2024, 1, 15, 0.0, "Jan"),
        (2024, 3, 20, 0.0, "Mar"),
        (2024, 6, 21, 0.0, "Jun"),
        (2024, 9, 22, 0.0, "Sep"),
        (2024, 12, 21, 0.0, "Dec"),
    ]

    latitudes = [0.0, 20.0, 40.0, 50.0, 55.0, 60.0, 63.0, 65.0]

    for y, m, d, h, dname in dates:
        jd = swe.julday(y, m, d, h)
        for lat in latitudes:
            geopos = (0.0, lat, 0.0)

            for rsmi, evname in [
                (CALC_RISE, "Rise"),
                (CALC_SET, "Set"),
            ]:
                label = f"Moon {evname} lat={lat}° {dname}"
                se_val = se_rise_trans(jd, MOON, geopos, rsmi)
                le_val = le_rise_trans(jd, MOON, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: both None")
                    continue
                if se_val is None or le_val is None:
                    r.fail(f"{label}: SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)

                # Moon: wider tolerance (faster motion, more sensitive to ephemeris)
                tol = 60.0 if lat <= 55.0 else 180.0
                if diff_s > tol:
                    r.fail(f"{label}: diff={diff_s:.1f}s")
                else:
                    r.ok(diff_s, label)

        print(f"  {dname}: tested {len(latitudes)} latitudes")

    return r.summary(), r


# ============================================================
# PART 4: Planets rise/set at extreme latitudes
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Planet rise/set at extreme latitudes")
    print("=" * 70)

    r = R("P4: Planet Rise/Set")

    jd = swe.julday(2024, 6, 21, 0.0)
    planets = [
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
    ]
    latitudes = [0.0, 30.0, 50.0, 60.0, 65.0]

    for lat in latitudes:
        geopos = (0.0, lat, 0.0)
        for body_id, body_name in planets:
            for rsmi, evname in [
                (CALC_RISE, "Rise"),
                (CALC_SET, "Set"),
            ]:
                label = f"{body_name} {evname} lat={lat}°"
                se_val = se_rise_trans(jd, body_id, geopos, rsmi)
                le_val = le_rise_trans(jd, body_id, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: both None")
                    continue
                if se_val is None or le_val is None:
                    r.fail(f"{label}: SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)
                tol = 60.0 if lat <= 50.0 else 180.0
                if diff_s > tol:
                    r.fail(f"{label}: diff={diff_s:.1f}s")
                else:
                    r.ok(diff_s, label)

        print(f"  lat={lat}°: tested {len(planets)} planets")

    return r.summary(), r


# ============================================================
# PART 5: Transit times across all latitudes
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Transit times (upper & lower) across latitudes")
    print("=" * 70)

    r = R("P5: Transits")

    jd = swe.julday(2024, 3, 20, 0.0)

    latitudes = [0.0, 20.0, 40.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0]
    bodies = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
    ]

    for lat in latitudes:
        geopos = (15.0, lat, 0.0)
        for body_id, body_name in bodies:
            for rsmi, evname in [
                (CALC_MTRANSIT, "UpperTr"),
                (CALC_ITRANSIT, "LowerTr"),
            ]:
                label = f"{body_name} {evname} lat={lat}°"
                se_val = se_rise_trans(jd, body_id, geopos, rsmi)
                le_val = le_rise_trans(jd, body_id, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: both None")
                    continue
                if se_val is None or le_val is None:
                    r.fail(f"{label}: SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)
                tol = 30.0  # Transits should always be precise
                if diff_s > tol:
                    r.fail(f"{label}: diff={diff_s:.1f}s")
                else:
                    r.ok(diff_s, label)

        print(f"  lat={lat}°: tested {len(bodies)} bodies")

    return r.summary(), r


# ============================================================
# PART 6: Twilight at extreme latitudes
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Twilight at extreme latitudes")
    print("=" * 70)

    r = R("P6: Twilight")

    jd = swe.julday(2024, 3, 20, 0.0)  # Equinox — twilight exists everywhere

    latitudes = [0.0, 30.0, 50.0, 60.0, 63.0, 65.0]

    twilight_flags = [
        (BIT_CIVIL_TWILIGHT, "Civil"),
        (BIT_NAUTIC_TWILIGHT, "Nautical"),
        (BIT_ASTRO_TWILIGHT, "Astronomical"),
    ]

    for lat in latitudes:
        geopos = (0.0, lat, 0.0)
        for tw_flag, tw_name in twilight_flags:
            for rsmi_base, evname in [
                (CALC_RISE, "Begin"),
                (CALC_SET, "End"),
            ]:
                rsmi = rsmi_base | tw_flag
                label = f"{tw_name} {evname} lat={lat}°"

                se_val = se_rise_trans(jd, SUN, geopos, rsmi)
                le_val = le_rise_trans(jd, SUN, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: both None (no twilight?)")
                    continue
                if se_val is None or le_val is None:
                    r.fail(f"{label}: SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)
                tol = 60.0 if lat <= 55 else 180.0
                if diff_s > tol:
                    r.fail(f"{label}: diff={diff_s:.1f}s")
                else:
                    r.ok(diff_s, label)

        print(f"  lat={lat}°: tested {len(twilight_flags)} twilight types")

    # Also test summer solstice — white nights at high lat
    print("  --- Summer solstice ---")
    jd_summer = swe.julday(2024, 6, 21, 0.0)
    for lat in [55.0, 58.0, 60.0, 63.0, 65.0]:
        geopos = (0.0, lat, 0.0)
        for tw_flag, tw_name in twilight_flags:
            for rsmi_base, evname in [
                (CALC_RISE, "Begin"),
                (CALC_SET, "End"),
            ]:
                rsmi = rsmi_base | tw_flag
                label = f"{tw_name} {evname} lat={lat}° Summer"

                se_val = se_rise_trans(jd_summer, SUN, geopos, rsmi)
                le_val = le_rise_trans(jd_summer, SUN, geopos, rsmi)

                if se_val is None and le_val is None:
                    r.skip(f"{label}: no twilight (white night)")
                    continue
                if se_val is None or le_val is None:
                    r.fail(f"{label}: DISAGREE SE={se_val} LE={le_val}")
                    continue

                diff_s = jd_diff_seconds(se_val, le_val)
                tol = 300.0  # Very wide for near-polar twilight
                if diff_s > tol:
                    r.fail(f"{label}: diff={diff_s:.1f}s")
                else:
                    r.ok(diff_s, label)

        print(f"  lat={lat}° Summer: tested")

    return r.summary(), r


# ============================================================
# PART 7: True horizon (elevated observer)
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: True horizon (elevated observer) at various latitudes")
    print("=" * 70)

    r = R("P7: True Horizon")

    jd = swe.julday(2024, 3, 20, 0.0)

    # Various elevations and latitudes
    test_cases = [
        (0.0, 40.0, 100.0, 0.0, "100m lat=40°"),
        (0.0, 40.0, 500.0, 0.0, "500m lat=40°"),
        (0.0, 40.0, 2000.0, 0.0, "2000m lat=40°"),
        (0.0, 60.0, 100.0, 0.0, "100m lat=60°"),
        (0.0, 60.0, 1000.0, 0.0, "1000m lat=60°"),
        # Mountain with a horizon altitude
        (0.0, 45.0, 1500.0, 2.0, "1500m horiz=2° lat=45°"),
        (0.0, 45.0, 0.0, -0.5, "sea horiz=-0.5° lat=45°"),
    ]

    for lon, lat, alt, horiz_alt, desc in test_cases:
        geopos = (lon, lat, alt)
        for rsmi, evname in [
            (CALC_RISE, "Rise"),
            (CALC_SET, "Set"),
        ]:
            label = f"Sun {evname} {desc}"
            se_val = se_rise_trans_true_hor(jd, SUN, geopos, rsmi, horiz_alt)
            le_val = le_rise_trans_true_hor(jd, SUN, geopos, rsmi, horiz_alt)

            if se_val is None and le_val is None:
                r.skip(f"{label}: both None")
                continue
            if se_val is None or le_val is None:
                r.fail(f"{label}: SE={se_val} LE={le_val}")
                continue

            diff_s = jd_diff_seconds(se_val, le_val)
            tol = 60.0
            if diff_s > tol:
                r.fail(f"{label}: diff={diff_s:.1f}s")
            else:
                r.ok(diff_s, label)

        print(f"  {desc}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 21: Rise/Set at Extreme Latitudes Deep Verification")
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
    print("ROUND 21 FINAL SUMMARY")
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
        print(f"\n  >>> ROUND 21: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 21: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
