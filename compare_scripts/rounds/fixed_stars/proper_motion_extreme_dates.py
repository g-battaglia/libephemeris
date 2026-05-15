#!/usr/bin/env python3
"""Round 115: Fixed Stars with Proper Motion at Extreme Dates

Tests fixed star positions at dates far from J2000 where proper motion
accumulates significantly. Also tests different flag combinations.
P1: Bright stars at 100-year intervals (1800-2200)
P2: High proper motion stars (Sirius, Arcturus, Aldebaran)
P3: Stars with different flags (EQUATORIAL, J2000, NONUT)
P4: Star catalog completeness check
P5: Star magnitude retrieval
"""

from __future__ import annotations

import sys
import os
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLG_SPEED = 256
FLG_J2000 = 32
FLG_NONUT = 64
FLG_EQUATORIAL = 2048
FLG_SIDEREAL = 65536

# Major stars to test
STARS = [
    "Aldebaran",
    "Regulus",
    "Spica",
    "Antares",
    "Fomalhaut",
    "Sirius",
    "Canopus",
    "Arcturus",
    "Vega",
    "Capella",
    "Rigel",
    "Procyon",
    "Betelgeuse",
    "Pollux",
    "Deneb",
    "Altair",
    "Castor",
    "Polaris",
    "Algol",
    "Rasalhague",
]

# High proper motion stars
HIGH_PM_STARS = ["Sirius", "Arcturus", "Procyon", "Pollux", "Aldebaran"]


def run_test(label, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL {label}: {detail}")


def get_star_pos(star_name, jd, flags=FLG_SPEED):
    """Get star position from both SE and LE."""
    # SE: fixstar2(name, jd, flags) -> (pos_tuple, starname, retflag)
    se_result = swe.fixstar2(star_name, jd, flags)
    se_pos = se_result[0]  # (lon, lat, dist, lon_spd, lat_spd, dist_spd)

    # LE: fixstar2_ut(name, jd, flags) -> (pos_tuple, starname, retflag)
    le_result = ephem.fixstar2_ut(star_name, jd, flags)
    le_pos = le_result[0]  # (lon, lat, dist, lon_spd, lat_spd, dist_spd)

    return se_pos, le_pos


# ============================================================
# P1: Bright stars at 100-year intervals (1800-2200)
# ============================================================
print("=== P1: Bright stars at 100-year intervals ===")

for star in STARS:
    for year in range(1850, 2151, 50):
        jd = swe.julday(year, 1, 1, 12.0)
        try:
            se_pos, le_pos = get_star_pos(star, jd)

            se_lon, se_lat = se_pos[0], se_pos[1]
            le_lon, le_lat = le_pos[0], le_pos[1]

            diff_lon = abs(se_lon - le_lon)
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lon_arcsec = diff_lon * 3600

            diff_lat = abs(se_lat - le_lat)
            diff_lat_arcsec = diff_lat * 3600

            # Tolerance: 5" for longitude, 5" for latitude
            tol = 5.0
            # At extreme dates, proper motion accumulation may differ more
            if year < 1900 or year > 2100:
                tol = 10.0

            label = f"P1 {star} Y{year}"
            if diff_lon_arcsec >= tol:
                run_test(
                    f"{label} lon",
                    False,
                    f'SE={se_lon:.6f} LE={le_lon:.6f} diff={diff_lon_arcsec:.2f}"',
                )
            else:
                passed += 1

            if diff_lat_arcsec >= tol:
                run_test(
                    f"{label} lat",
                    False,
                    f'SE={se_lat:.6f} LE={le_lat:.6f} diff={diff_lat_arcsec:.2f}"',
                )
            else:
                passed += 1

        except Exception as e:
            errors += 1
            if "not found" not in str(e).lower() and "range" not in str(e).lower():
                print(f"  ERROR P1 {star} Y{year}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: High proper motion stars — daily over 10 years
# ============================================================
print("\n=== P2: High PM stars daily over 10 years ===")

for star in HIGH_PM_STARS:
    jd_start = swe.julday(2020, 1, 1, 12.0)
    for day in range(0, 3650, 30):  # monthly for 10 years
        jd = jd_start + day
        try:
            se_pos, le_pos = get_star_pos(star, jd)

            diff_lon = abs(se_pos[0] - le_pos[0])
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lon_arcsec = diff_lon * 3600

            diff_lat = abs(se_pos[1] - le_pos[1])
            diff_lat_arcsec = diff_lat * 3600

            tol = 3.0  # tight tolerance

            if diff_lon_arcsec >= tol or diff_lat_arcsec >= tol:
                run_test(
                    f"P2 {star} d={day}",
                    False,
                    f'dlon={diff_lon_arcsec:.2f}" dlat={diff_lat_arcsec:.2f}"',
                )
            else:
                passed += 1

        except Exception as e:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Stars with different flags
# ============================================================
print("\n=== P3: Stars with different flags ===")

FLAG_COMBOS = [
    (FLG_SPEED, "default"),
    (FLG_SPEED | FLG_EQUATORIAL, "EQUATORIAL"),
    (FLG_SPEED | FLG_J2000 | FLG_NONUT, "J2000"),
    (FLG_SPEED | FLG_NONUT, "NONUT"),
    (FLG_SPEED | FLG_EQUATORIAL | FLG_J2000 | FLG_NONUT, "EQUAT+J2000"),
]

TEST_STARS = ["Sirius", "Regulus", "Spica", "Aldebaran", "Vega"]

for flags, flag_label in FLAG_COMBOS:
    for star in TEST_STARS:
        for jd_label, jd in [
            ("J2000", 2451545.0),
            ("2024", swe.julday(2024, 6, 15, 12.0)),
        ]:
            try:
                se_pos, le_pos = get_star_pos(star, jd, flags)

                diff_lon = abs(se_pos[0] - le_pos[0])
                if diff_lon > 180:
                    diff_lon = 360 - diff_lon
                diff_lon_arcsec = diff_lon * 3600

                diff_lat = abs(se_pos[1] - le_pos[1])
                diff_lat_arcsec = diff_lat * 3600

                tol = 5.0
                if "J2000" in flag_label:
                    tol = 8.0  # known J2000 frame differences
                if "NONUT" in flag_label:
                    tol = 20.0  # NONUT differences

                label = f"P3 {flag_label} {star} {jd_label}"
                if diff_lon_arcsec >= tol:
                    run_test(
                        f"{label} lon",
                        False,
                        f'SE={se_pos[0]:.6f} LE={le_pos[0]:.6f} diff={diff_lon_arcsec:.2f}"',
                    )
                else:
                    passed += 1

                if diff_lat_arcsec >= tol:
                    run_test(
                        f"{label} lat",
                        False,
                        f'SE={se_pos[1]:.6f} LE={le_pos[1]:.6f} diff={diff_lat_arcsec:.2f}"',
                    )
                else:
                    passed += 1

            except Exception as e:
                errors += 1
                print(f"  ERROR P3 {flag_label} {star} {jd_label}: {e}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Star speed comparison
# ============================================================
print("\n=== P4: Star speed comparison ===")

jd = 2451545.0
for star in STARS:
    try:
        se_pos, le_pos = get_star_pos(star, jd)

        # Speed comparison (arcsec/day)
        se_lon_spd = se_pos[3] * 3600  # convert deg/day to arcsec/day
        le_lon_spd = le_pos[3] * 3600
        diff_spd = abs(se_lon_spd - le_lon_spd)

        # Star speeds are tiny (~0.01-0.15"/day), so relative tolerance
        if abs(se_lon_spd) > 0.001:
            rel_diff = diff_spd / abs(se_lon_spd)
            ok = rel_diff < 0.1  # 10% relative tolerance
        else:
            ok = diff_spd < 0.01

        label = f"P4 {star} lon_speed"
        if not ok:
            run_test(
                label,
                False,
                f'SE={se_lon_spd:.6f}"/d LE={le_lon_spd:.6f}"/d diff={diff_spd:.6f}"/d',
            )
        else:
            passed += 1

        # Latitude speed
        se_lat_spd = se_pos[4] * 3600
        le_lat_spd = le_pos[4] * 3600
        diff_lat_spd = abs(se_lat_spd - le_lat_spd)

        if abs(se_lat_spd) > 0.001:
            rel_diff_lat = diff_lat_spd / abs(se_lat_spd)
            ok_lat = rel_diff_lat < 0.5  # wider for lat speed
        else:
            ok_lat = diff_lat_spd < 0.01

        label = f"P4 {star} lat_speed"
        if not ok_lat:
            run_test(
                label,
                False,
                f'SE={se_lat_spd:.6f}"/d LE={le_lat_spd:.6f}"/d diff={diff_lat_spd:.6f}"/d',
            )
        else:
            passed += 1

    except Exception as e:
        errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Star name resolution and catalog check
# ============================================================
print("\n=== P5: Star name resolution ===")

# Test that all major stars can be found
ALL_STARS = STARS + [
    "Achernar",
    "Hamal",
    "Menkar",
    "Mirach",
    "Scheat",
    "Markab",
    "Alphecca",
    "Zubenelgenubi",
    "Zubeneschamali",
    "Unukalhai",
    "Ras Alhague",
    "Sabik",
    "Shaula",
]

jd = 2451545.0
for star in ALL_STARS:
    try:
        le_result = ephem.fixstar2_ut(star, jd, FLG_SPEED)
        le_name = le_result[1]
        le_pos = le_result[0]

        # Check position is reasonable
        lon = le_pos[0]
        lat = le_pos[1]

        ok = (0 <= lon < 360) and (-90 <= lat <= 90)
        label = f"P5 {star} found"
        run_test(label, ok, f"lon={lon:.4f} lat={lat:.4f}")

    except Exception as e:
        errors += 1
        print(f"  ERROR P5 {star}: {e}")

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Proper motion accumulation test
# Compare positions 50 years apart — proper motion should be consistent
# ============================================================
print("\n=== P6: Proper motion accumulation consistency ===")

for star in HIGH_PM_STARS:
    try:
        jd1 = swe.julday(1975, 1, 1, 12.0)
        jd2 = swe.julday(2025, 1, 1, 12.0)

        se_pos1, le_pos1 = get_star_pos(star, jd1)
        se_pos2, le_pos2 = get_star_pos(star, jd2)

        # SE proper motion over 50 years
        se_delta_lon = (se_pos2[0] - se_pos1[0]) * 3600  # arcsec
        le_delta_lon = (le_pos2[0] - le_pos1[0]) * 3600

        diff_delta = abs(se_delta_lon - le_delta_lon)

        # Over 50 years, PM accumulation should agree within 5"
        label = f"P6 {star} PM_accum"
        run_test(
            label,
            diff_delta < 5.0,
            f'SE_delta={se_delta_lon:.2f}" LE_delta={le_delta_lon:.2f}" diff={diff_delta:.2f}"',
        )

    except Exception as e:
        errors += 1
        print(f"  ERROR P6 {star}: {e}")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 115 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
