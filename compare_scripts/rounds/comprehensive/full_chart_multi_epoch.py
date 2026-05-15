#!/usr/bin/env python3
"""Round 77: Full Chart Integration Multi-Epoch"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256  # FLG_SPEED


def se_hsys(ch):
    return ch.encode("ascii")


def le_hsys(ch):
    return ord(ch)


# Famous chart dates (approximate JDs)
CHARTS = [
    (2451545.0, 41.90, 12.50, "J2000 Rome"),
    (2415020.0, 48.85, 2.35, "1900 Paris"),
    (2440587.5, 40.78, -73.97, "1970 NYC"),
    (2460000.0, 35.69, 139.69, "2023 Tokyo"),
    (2430000.0, 51.51, -0.13, "1941 London"),
    (2445000.0, -33.87, 151.21, "1982 Sydney"),
    (2455000.0, 28.61, 77.21, "2009 Delhi"),
    (2437000.0, 55.75, 37.62, "1960 Moscow"),
    (2450000.0, 19.43, -99.13, "1995 Mexico City"),
    (2442000.0, -22.91, -43.17, "1974 Rio"),
    (2448000.0, 59.33, 18.07, "1990 Stockholm"),
    (2453000.0, 1.35, 103.82, "2004 Singapore"),
]

HOUSE_SYSTEMS = ["P", "K", "O", "R", "C", "E", "W", "B"]

print("=" * 70)
print("ROUND 77: Full Chart Integration Multi-Epoch")
print("=" * 70)

# P1: Complete chart — all planets + houses
print("\n=== P1: Full chart for all epochs ===")
for jd, lat, lon, name in CHARTS:
    chart_pass = chart_fail = 0
    # Planets
    for body in range(10):
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 1.0:
                passed += 1
                chart_pass += 1
            else:
                failed += 1
                chart_fail += 1
        except:
            errors += 1
    # Nodes
    for body in [10, 11, 12, 13]:
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 1.0:
                passed += 1
                chart_pass += 1
            else:
                failed += 1
                chart_fail += 1
        except:
            errors += 1
    # Houses
    for hsys in HOUSE_SYSTEMS:
        try:
            se_h = swe.houses_ex(jd, lat, lon, se_hsys(hsys))
            le_h = ephem.houses_ex(jd, lat, lon, le_hsys(hsys), 0)
            for i in range(12):
                diff = abs(se_h[0][i] - le_h[0][i])
                if diff > 180:
                    diff = 360 - diff
                if diff * 3600 < 5.0:
                    passed += 1
                    chart_pass += 1
                else:
                    failed += 1
                    chart_fail += 1
        except:
            errors += 1
    if chart_fail > 0:
        print(f"  {name}: {chart_pass}/{chart_pass + chart_fail}")
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Ascendant/MC consistency
print("\n=== P2: Asc/MC consistency ===")
for jd, lat, lon, name in CHARTS:
    try:
        se_h = swe.houses_ex(jd, lat, lon, se_hsys("P"))
        le_h = ephem.houses_ex(jd, lat, lon, le_hsys("P"), 0)
        # ASC
        diff_asc = abs(se_h[1][0] - le_h[1][0])
        if diff_asc > 180:
            diff_asc = 360 - diff_asc
        if diff_asc * 3600 < 1.0:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL ASC {name}: SE={se_h[1][0]:.4f} LE={le_h[1][0]:.4f}")
        # MC
        diff_mc = abs(se_h[1][1] - le_h[1][1])
        if diff_mc > 180:
            diff_mc = 360 - diff_mc
        if diff_mc * 3600 < 1.0:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL MC {name}: SE={se_h[1][1]:.4f} LE={le_h[1][1]:.4f}")
    except:
        errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Sidereal time
print("\n=== P3: Sidereal time ===")
for jd, lat, lon, name in CHARTS:
    try:
        se_st = swe.sidtime(jd)
        le_st = ephem.sidtime(jd)
        diff = abs(se_st - le_st) * 3600  # hours to seconds
        if diff < 0.1:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL ST {name}: SE={se_st:.8f} LE={le_st:.8f} diff={diff:.4f}s")
    except:
        errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Delta-T consistency
print("\n=== P4: Delta-T ===")
for jd, lat, lon, name in CHARTS:
    try:
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff = abs(se_dt - le_dt) * 86400  # days to seconds
        if diff < 5.0:  # within 5 seconds
            passed += 1
        else:
            failed += 1
            print(f"  FAIL DT {name}: SE={se_dt:.10f} LE={le_dt:.10f} diff={diff:.2f}s")
    except:
        errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: All planets + all 6 position values at J2000
print("\n=== P5: Full 6-value comparison at J2000 ===")
jd = 2451545.0
for body in range(10):
    try:
        se = swe.calc_ut(jd, body, FLAGS)
        le = ephem.calc_ut(jd, body, FLAGS)
        for idx in range(6):
            diff = abs(se[0][idx] - le[0][idx])
            if idx == 0 and diff > 180:
                diff = 360 - diff
            if idx < 3:
                tol = 0.001 if idx < 2 else 0.0001
            else:
                tol = 0.001
            if diff < tol:
                passed += 1
            else:
                failed += 1
    except:
        errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 77 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
