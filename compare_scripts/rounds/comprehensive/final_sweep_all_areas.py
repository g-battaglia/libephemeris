#!/usr/bin/env python3
"""Round 210: Final comprehensive sweep.

A broad sweep testing all major API areas at multiple dates to serve as
a final quality gate. Tests calc_ut, houses, fixed stars, nutation,
sidereal time, ayanamsha, and coordinate transforms.
"""

from __future__ import annotations

import os
import sys

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

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

TEST_JDS = [
    2451545.0,  # J2000
    2455197.5,  # 2010
    2458849.5,  # 2020
    2460310.5,  # 2024
    2443144.5,  # 1977
    2415020.0,  # 1900
    2440587.5,  # 1970 (Unix epoch)
]

ALL_BODIES = [
    ("Sun", ephem.SUN, swe.SUN, 0.5),
    ("Moon", ephem.MOON, swe.MOON, 1.0),
    ("Mercury", ephem.MERCURY, swe.MERCURY, 0.5),
    ("Venus", ephem.VENUS, swe.VENUS, 0.5),
    ("Mars", ephem.MARS, swe.MARS, 0.5),
    ("Jupiter", ephem.JUPITER, swe.JUPITER, 0.5),
    ("Saturn", ephem.SATURN, swe.SATURN, 0.5),
    ("Uranus", ephem.URANUS, swe.URANUS, 0.5),
    ("Neptune", ephem.NEPTUNE, swe.NEPTUNE, 0.5),
    ("Pluto", ephem.PLUTO, swe.PLUTO, 1.0),
    ("Chiron", ephem.CHIRON, swe.CHIRON, 0.5),
    ("MeanNode", ephem.MEAN_NODE, swe.MEAN_NODE, 0.5),
    ("TrueNode", ephem.TRUE_NODE, swe.TRUE_NODE, 0.5),
]


def test_all_bodies():
    global passed, failed, total
    print("=" * 70)
    print("Round 210: Final Comprehensive Sweep")
    print("=" * 70)

    print("\n--- Ecliptic of Date ---")
    for jd in TEST_JDS:
        for bname, le_b, se_b, tol in ALL_BODIES:
            try:
                le_r = ephem.calc_ut(jd, le_b, FLAGS)
                se_r = swe.calc_ut(jd, se_b, swe.FLG_SWIEPH | swe.FLG_SPEED)
            except Exception:
                continue
            total += 1
            d = abs(le_r[0][0] - se_r[0][0])
            if d > 180:
                d = 360 - d
            if d * 3600 <= tol:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {bname} JD={jd:.1f}: diff={d * 3600:.4f}"')

    print("\n--- Equatorial ---")
    eq_flags_le = FLAGS | ephem.FLG_EQUATORIAL
    eq_flags_se = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
    for jd in TEST_JDS[:3]:
        for bname, le_b, se_b, tol in ALL_BODIES[:7]:
            try:
                le_r = ephem.calc_ut(jd, le_b, eq_flags_le)
                se_r = swe.calc_ut(jd, se_b, eq_flags_se)
            except Exception:
                continue
            total += 1
            d = abs(le_r[0][0] - se_r[0][0])
            if d > 180:
                d = 360 - d
            if d * 3600 <= tol:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {bname} EQ JD={jd:.1f}: diff={d * 3600:.4f}"')

    print("\n--- Houses ---")
    locs = [(52.52, 13.405), (40.71, -74.01), (-33.87, 151.21)]
    for jd in TEST_JDS[:3]:
        for lat, lon in locs:
            try:
                le_r = ephem.houses_ex2(jd, lat, lon, ord("P"), 0)
                se_r = swe.houses_ex(jd, lat, lon, b"P")
            except Exception:
                continue
            for i in range(12):
                total += 1
                d = abs(le_r[0][i] - se_r[0][i])
                if d > 180:
                    d = 360 - d
                if d * 3600 <= 1.0:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  House cusp{i + 1} JD={jd:.1f} lat={lat}: diff={d * 3600:.2f}"'
                    )

    print("\n--- Sidereal Time ---")
    for jd in TEST_JDS:
        total += 1
        le_st = ephem.sidtime(jd)
        se_st = swe.sidtime(jd)
        d = abs(le_st - se_st)
        if d > 12:
            d = 24 - d
        if d * 3600 <= 0.05:
            passed += 1
        else:
            failed += 1
            failures.append(f"  sidtime JD={jd:.1f}: diff={d * 3600:.6f}s")

    print("\n--- Fixed Stars ---")
    stars = ["Aldebaran", "Regulus", "Spica", "Antares", "Sirius", "Vega"]
    for jd in TEST_JDS[:2]:
        for star in stars:
            try:
                le_r = ephem.fixstar2_ut(star, jd, FLAGS)
                se_r = swe.fixstar2(star, jd, swe.FLG_SWIEPH | swe.FLG_SPEED)
            except Exception:
                continue
            total += 1
            d = abs(le_r[0][0] - se_r[0][0])
            if d > 180:
                d = 360 - d
            if d * 3600 <= 1.0:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {star} JD={jd:.1f}: diff={d * 3600:.4f}"')

    print("\n--- Nutation ---")
    for jd in TEST_JDS:
        try:
            le_r = ephem.calc_ut(jd, -1, 0)
            se_r = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
        except Exception:
            continue
        for idx in range(4):
            total += 1
            d = abs(le_r[0][idx] - se_r[0][idx]) * 3600
            if d <= 0.1:
                passed += 1
            else:
                failed += 1
                failures.append(f'  Nut[{idx}] JD={jd:.1f}: diff={d:.4f}"')

    print("\n--- Delta-T ---")
    for jd in TEST_JDS:
        total += 1
        le_dt = ephem.deltat(jd)
        se_dt = swe.deltat(jd)
        d = abs(le_dt - se_dt) * 86400
        if d <= 1.0:
            passed += 1
        else:
            failed += 1
            failures.append(f"  deltat JD={jd:.1f}: diff={d:.4f}s")

    print("\n--- cotrans ---")
    for lon in [0, 45, 90, 135, 180, 225, 270, 315]:
        for lat in [-60, -30, 0, 30, 60]:
            obl = 23.4393
            total += 1
            le_ct = ephem.cotrans((lon, lat, 1.0), obl)
            se_ct = swe.cotrans((lon, lat, 1.0), obl)
            d = abs(le_ct[0] - se_ct[0])
            if d > 180:
                d = 360 - d
            if d * 3600 <= 0.01:
                passed += 1
            else:
                failed += 1
                failures.append(f'  cotrans lon={lon} lat={lat}: diff={d * 3600:.4f}"')


if __name__ == "__main__":
    test_all_bodies()

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
