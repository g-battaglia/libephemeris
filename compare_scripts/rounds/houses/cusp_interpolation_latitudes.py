#!/usr/bin/env python3
"""Round 110: House Cusp Interpolation Across Latitudes

Deep test of house cusps across a fine latitude grid for multiple house systems.
P1: Placidus cusps at 5° latitude intervals (-60 to +60)
P2: Koch cusps at 5° latitude intervals
P3: Regiomontanus cusps at 5° latitude intervals
P4: Campanus cusps at 5° latitude intervals
P5: Equal/Whole Sign cusps (should be latitude-independent for Equal)
P6: Cusps at tropical latitudes (23.44° N/S) — seasonal extremes
P7: Cusp continuity — small latitude steps, cusps should change smoothly
P8: Multiple ARMC values at fixed latitude
P9: Porphyry cusps deep
P10: Morinus cusps (should be latitude-independent)
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLG_SPEED = 256


def se_hsys(ch):
    return ch.encode("ascii")


def le_hsys(ch):
    return ord(ch)


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


def compare_cusps(se_cusps, le_cusps, label, tol_arcsec=5.0, ncusps=12):
    """Compare house cusps, return (n_pass, n_fail)."""
    global passed, failed
    np, nf = 0, 0
    for i in range(ncusps):
        se_c = se_cusps[i]
        le_c = le_cusps[i]
        diff = abs(se_c - le_c)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600

        if diff_arcsec < tol_arcsec:
            passed += 1
            np += 1
        else:
            failed += 1
            nf += 1
            print(
                f'  FAIL {label} cusp{i}: SE={se_c:.6f} LE={le_c:.6f} diff={diff_arcsec:.2f}"'
            )
    return np, nf


def compare_ascmc(se_ascmc, le_ascmc, label, tol_arcsec=5.0, indices=None):
    """Compare ascmc values."""
    global passed, failed
    if indices is None:
        indices = [0, 1, 2, 3]  # ASC, MC, ARMC, Vertex
    for i in indices:
        names = {0: "ASC", 1: "MC", 2: "ARMC", 3: "Vertex"}
        name = names.get(i, f"idx{i}")
        se_v = se_ascmc[i]
        le_v = le_ascmc[i]
        diff = abs(se_v - le_v)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600

        if diff_arcsec < tol_arcsec:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL {label} {name}: SE={se_v:.6f} LE={le_v:.6f} diff={diff_arcsec:.2f}"'
            )


# Test epochs
TEST_JDS = [
    (2451545.0, "J2000"),
    (swe.julday(2024, 3, 20, 12.0), "2024-equi"),
    (swe.julday(2024, 6, 21, 12.0), "2024-sols"),
    (swe.julday(2024, 9, 22, 12.0), "2024-equ2"),
]


# ============================================================
# P1: Placidus cusps at 5° latitude intervals
# ============================================================
print("=== P1: Placidus cusps across latitudes ===")

for jd, epoch_label in TEST_JDS:
    for lat in range(-60, 61, 5):
        lon = 0.0
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("P"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("P"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(se_cusps, le_cusps, f"P1 Placidus lat={lat} {epoch_label}")
            compare_ascmc(se_ascmc, le_ascmc, f"P1 Placidus lat={lat} {epoch_label}")
        except Exception as e:
            errors += 1
            if "PolarCircle" not in str(e) and "polar" not in str(e).lower():
                print(f"  ERROR P1 lat={lat} {epoch_label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Koch cusps at 5° latitude intervals
# ============================================================
print("\n=== P2: Koch cusps across latitudes ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-55, 56, 5):
        lon = 11.25  # Rome-ish
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("K"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("K"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(se_cusps, le_cusps, f"P2 Koch lat={lat} {epoch_label}")
        except Exception as e:
            errors += 1
            if "PolarCircle" not in str(e) and "polar" not in str(e).lower():
                print(f"  ERROR P2 lat={lat} {epoch_label}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Regiomontanus cusps
# ============================================================
print("\n=== P3: Regiomontanus cusps across latitudes ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-65, 66, 5):
        lon = -73.97  # NYC
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("R"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("R"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(se_cusps, le_cusps, f"P3 Regio lat={lat} {epoch_label}")
        except Exception as e:
            errors += 1
            if "PolarCircle" not in str(e) and "polar" not in str(e).lower():
                print(f"  ERROR P3 lat={lat} {epoch_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Campanus cusps
# ============================================================
print("\n=== P4: Campanus cusps across latitudes ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-65, 66, 5):
        lon = 139.69  # Tokyo
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("C"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("C"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(se_cusps, le_cusps, f"P4 Campanus lat={lat} {epoch_label}")
        except Exception as e:
            errors += 1
            print(f"  ERROR P4 lat={lat} {epoch_label}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Equal house cusps (should be same at all latitudes)
# ============================================================
print("\n=== P5: Equal house cusps ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-80, 81, 10):
        lon = 0.0
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("E"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("E"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(se_cusps, le_cusps, f"P5 Equal lat={lat} {epoch_label}")
        except Exception as e:
            errors += 1
            print(f"  ERROR P5 lat={lat} {epoch_label}: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Cusps at tropical latitudes (23.44° N/S)
# ============================================================
print("\n=== P6: Cusps at tropical latitudes ===")

TROPICAL_LATS = [23.44, -23.44, 0.0, 45.0, -45.0]
HSYS_LIST = [("P", "Placidus"), ("K", "Koch"), ("R", "Regio"), ("C", "Campanus")]

for jd, epoch_label in TEST_JDS:
    for lat in TROPICAL_LATS:
        for hsys_ch, hsys_name in HSYS_LIST:
            try:
                se_result = swe.houses_ex(jd, lat, 0.0, se_hsys(hsys_ch))
                le_result = ephem.houses_ex(jd, lat, 0.0, le_hsys(hsys_ch), 0)

                se_cusps, se_ascmc = se_result[0], se_result[1]
                le_cusps, le_ascmc = le_result[0], le_result[1]

                compare_cusps(
                    se_cusps,
                    le_cusps,
                    f"P6 {hsys_name} lat={lat} {epoch_label}",
                )
            except Exception as e:
                errors += 1
                if "PolarCircle" not in str(e):
                    print(f"  ERROR P6 {hsys_name} lat={lat} {epoch_label}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Cusp continuity — small latitude steps
# ============================================================
print("\n=== P7: Cusp continuity (1° steps around 45°N) ===")

jd = 2451545.0  # J2000
prev_cusps = None

for lat_10 in range(350, 551):  # 35.0 to 55.0 in 0.1° steps
    lat = lat_10 / 10.0
    try:
        le_result = ephem.houses_ex(jd, lat, 0.0, le_hsys("P"), 0)
        le_cusps = le_result[0]

        if prev_cusps is not None:
            # Each cusp should change smoothly (up to ~5° per 0.1° lat step is normal)
            for i in range(12):
                diff = abs(le_cusps[i] - prev_cusps[i])
                if diff > 180:
                    diff = 360 - diff
                if diff > 5.0:
                    run_test(
                        f"P7 Placidus cusp{i} lat={lat:.1f}",
                        False,
                        f"jump={diff:.4f}° from lat={lat - 0.1:.1f}",
                    )
                else:
                    passed += 1

        prev_cusps = le_cusps
    except Exception as e:
        errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Multiple ARMC values at fixed latitude
# ============================================================
print("\n=== P8: Multiple ARMC values (houses_armc) ===")

obliquity = 23.4393  # Approximate for J2000

for lat in [0, 30, 45, 52, -33]:
    for armc in range(0, 360, 15):
        for hsys_ch, hsys_name in [
            ("P", "Placidus"),
            ("R", "Regio"),
            ("C", "Campanus"),
        ]:
            try:
                se_result = swe.houses_armc(
                    float(armc), float(lat), obliquity, se_hsys(hsys_ch)
                )
                le_result = ephem.houses_armc(
                    float(armc), float(lat), obliquity, le_hsys(hsys_ch)
                )

                se_cusps, se_ascmc = se_result[0], se_result[1]
                le_cusps, le_ascmc = le_result[0], le_result[1]

                compare_cusps(
                    se_cusps,
                    le_cusps,
                    f"P8 {hsys_name} armc={armc} lat={lat}",
                    tol_arcsec=2.0,
                )
            except Exception as e:
                errors += 1
                if "PolarCircle" not in str(e) and "polar" not in str(e).lower():
                    print(f"  ERROR P8 {hsys_name} armc={armc} lat={lat}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P9: Porphyry cusps
# ============================================================
print("\n=== P9: Porphyry cusps ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-70, 71, 10):
        lon = 0.0
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("O"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("O"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(
                se_cusps,
                le_cusps,
                f"P9 Porphyry lat={lat} {epoch_label}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P9 lat={lat} {epoch_label}: {e}")

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P10: Morinus cusps (latitude-independent)
# ============================================================
print("\n=== P10: Morinus cusps ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-80, 81, 20):
        lon = 0.0
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("M"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("M"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(
                se_cusps,
                le_cusps,
                f"P10 Morinus lat={lat} {epoch_label}",
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P10 lat={lat} {epoch_label}: {e}")

print(f"  After P10: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P11: Whole Sign houses
# ============================================================
print("\n=== P11: Whole Sign cusps ===")

for jd, epoch_label in TEST_JDS[:2]:
    for lat in range(-60, 61, 15):
        lon = 0.0
        try:
            se_result = swe.houses_ex(jd, lat, lon, se_hsys("W"))
            le_result = ephem.houses_ex(jd, lat, lon, le_hsys("W"), 0)

            se_cusps, se_ascmc = se_result[0], se_result[1]
            le_cusps, le_ascmc = le_result[0], le_result[1]

            compare_cusps(
                se_cusps,
                le_cusps,
                f"P11 WholeSgn lat={lat} {epoch_label}",
                tol_arcsec=1.0,
            )
        except Exception as e:
            errors += 1
            print(f"  ERROR P11 lat={lat} {epoch_label}: {e}")

print(f"  After P11: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 110 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
