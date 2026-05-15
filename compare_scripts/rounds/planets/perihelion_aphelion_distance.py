#!/usr/bin/env python3
"""Round 224: Planetary distance at perihelion/aphelion.

Finds heliocentric distance extremes for all planets and compares
positions at those moments.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

LE_FLAGS_H = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_HELCTR
SE_FLAGS_H = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_HELCTR
LE_FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED
SE_FLAGS = swe.FLG_SWIEPH | swe.FLG_SPEED

BODIES = [
    ("Mercury", ephem.MERCURY, swe.MERCURY, 10.0, 20),
    ("Venus", ephem.VENUS, swe.VENUS, 20.0, 12),
    ("Mars", ephem.MARS, swe.MARS, 30.0, 8),
    ("Jupiter", ephem.JUPITER, swe.JUPITER, 100.0, 6),
    ("Saturn", ephem.SATURN, swe.SATURN, 200.0, 4),
    ("Earth", ephem.EARTH, swe.EARTH, 5.0, 10),
]

JD_START = 2451545.0
JD_END = JD_START + 365.25 * 25


def find_dist_extremes(se_b, step, max_count):
    extremes = []
    jd = JD_START
    prev, pprev = None, None
    while jd < JD_END and len(extremes) < max_count:
        try:
            r = swe.calc_ut(jd, se_b, SE_FLAGS_H)
            d = r[0][2]
        except:
            jd += step
            continue
        if prev is not None and pprev is not None:
            if prev > pprev and prev > d:
                extremes.append(("Aphelion", jd - step))
            elif prev < pprev and prev < d:
                extremes.append(("Perihelion", jd - step))
        pprev, prev = prev, d
        jd += step
    return extremes


def compare_at(label, le_b, se_b, jd):
    global passed, failed, total
    for flags_label, le_f, se_f in [
        ("Helio", LE_FLAGS_H, SE_FLAGS_H),
        ("Geo", LE_FLAGS, SE_FLAGS),
    ]:
        try:
            le_r = ephem.calc_ut(jd, le_b, le_f)
            se_r = swe.calc_ut(jd, se_b, se_f)
        except:
            continue

        total += 1
        d = abs(le_r[0][0] - se_r[0][0])
        if d > 180:
            d = 360 - d
        tol = 2.0
        if d * 3600 <= tol:
            passed += 1
        else:
            failed += 1
            failures.append(f'  {label} {flags_label} LON: diff={d * 3600:.4f}"')

        total += 1
        dd = abs(le_r[0][2] - se_r[0][2])
        if dd <= 0.0001:
            passed += 1
        else:
            failed += 1
            failures.append(f"  {label} {flags_label} DIST: diff={dd:.8f} AU")


if __name__ == "__main__":
    print("=" * 70)
    print("Round 224: Planetary Distance at Perihelion/Aphelion")
    print("=" * 70)
    for bname, le_b, se_b, step, cnt in BODIES:
        print(f"\n--- {bname} ---")
        extremes = find_dist_extremes(se_b, step, cnt)
        print(f"  Found {len(extremes)} extremes")
        for etype, jd in extremes:
            compare_at(f"{bname} {etype}", le_b, se_b, jd)

    print(f"\n{'=' * 70}")
    if total > 0:
        print(
            f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
        )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
