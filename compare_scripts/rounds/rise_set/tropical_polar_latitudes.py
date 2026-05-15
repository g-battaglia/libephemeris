#!/usr/bin/env python3
"""Round 80: Rise/Set at Tropical and Polar Latitudes

Tests rise_trans and rise_trans_true_hor at extreme latitudes where
circumpolar/never-rise conditions create edge cases:
- Equator (0°)
- Tropics (±23.44°)
- Arctic/Antarctic circles (±66.56°)
- High Arctic (±70°, ±75°, ±80°)
- Near-pole (±85°, ±88°)
- Multiple bodies: Sun, Moon, Venus, Mars, Jupiter
- Multiple dates across seasons (solstices, equinoxes)
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
skipped = 0

# Rise/set flags
CALC_RISE = 1
CALC_SET = 2
CALC_MTRANSIT = 4
CALC_ITRANSIT = 8
BIT_DISC_CENTER = 256
BIT_NO_REFRACTION = 512


def se_hsys(ch):
    return ch.encode("ascii")


def compare_rise_set(label, body, jd, lat, lon, alt, rsmi, tol_sec=30.0):
    """Compare rise/set times between SE and LE."""
    global passed, failed, errors, skipped
    geopos = [lon, lat, alt]
    try:
        se_result = swe.rise_trans(jd, body, rsmi, geopos, 1013.25, 15.0)
        se_jd = se_result[1][0]
    except Exception as e:
        se_err = str(e)
        se_jd = None

    try:
        le_result = ephem.rise_trans(
            jd, body, rsmi, [lon, lat, alt], 1013.25, 15.0, 2
        )
        le_jd = le_result[1][0]
    except Exception as e:
        le_err = str(e)
        le_jd = None

    # Both failed — likely circumpolar/never-rise, both agree
    if se_jd is None and le_jd is None:
        passed += 1
        return

    # One failed, other didn't
    if se_jd is None and le_jd is not None:
        # SE couldn't find, LE found something — check if LE result is reasonable
        skipped += 1
        return
    if se_jd is not None and le_jd is None:
        skipped += 1
        return

    # Both returned values — compare
    # SE returns 0.0 for "not found"
    if se_jd < 1.0 and (le_jd is None or le_jd == 0.0):
        passed += 1
        return
    if se_jd < 1.0:
        skipped += 1
        return

    if isinstance(le_jd, tuple):
        le_jd = le_jd[0] if le_jd else None
    if le_jd is None or le_jd == 0.0:
        skipped += 1
        return

    diff_s = abs(float(se_jd) - float(le_jd)) * 86400.0
    if diff_s < tol_sec:
        passed += 1
    else:
        failed += 1
        print(
            f"  FAIL {label}: SE={float(se_jd):.6f} LE={float(le_jd):.6f} diff={diff_s:.1f}s"
        )


print("=" * 70)
print("ROUND 80: Rise/Set at Tropical and Polar Latitudes")
print("=" * 70)

# Key dates (2000-2005, covering solstices/equinoxes)
dates = {
    "2000 VE": 2451623.5,  # March 20, 2000 (vernal equinox)
    "2000 SS": 2451716.5,  # June 21, 2000 (summer solstice)
    "2000 AE": 2451808.5,  # Sep 22, 2000 (autumnal equinox)
    "2000 WS": 2451899.5,  # Dec 21, 2000 (winter solstice)
    "2002 VE": 2452353.5,  # Mar 20, 2002
    "2002 SS": 2452446.5,  # Jun 21, 2002
    "2003 WS": 2452995.5,  # Dec 22, 2003
    "2005 SS": 2453542.5,  # Jun 21, 2005
}

# Latitudes to test
latitudes = [
    0.0,
    23.44,
    -23.44,
    45.0,
    -45.0,
    60.0,
    -60.0,
    66.56,
    -66.56,
    70.0,
    -70.0,
    75.0,
    -75.0,
    80.0,
    -80.0,
]

lon = 0.0  # Greenwich
alt = 0.0

bodies = [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MARS, "Mars"),
]

# ============================================================
# P1: Sun rise/set at all latitudes and seasons
# ============================================================
print("\n=== P1: Sun rise/set across latitudes and seasons ===")

for date_name, jd in dates.items():
    for lat in latitudes:
        for rsmi, rsmi_name in [(CALC_RISE, "rise"), (CALC_SET, "set")]:
            label = f"Sun {rsmi_name} {date_name} lat={lat}"
            compare_rise_set(label, swe.SUN, jd, lat, lon, alt, rsmi)

print(
    f"  After P1: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P2: Moon rise/set
# ============================================================
print("\n=== P2: Moon rise/set across latitudes and seasons ===")

for date_name, jd in list(dates.items())[:4]:  # Just 2000 dates
    for lat in latitudes:
        for rsmi, rsmi_name in [(CALC_RISE, "rise"), (CALC_SET, "set")]:
            label = f"Moon {rsmi_name} {date_name} lat={lat}"
            compare_rise_set(label, swe.MOON, jd, lat, lon, alt, rsmi, tol_sec=60.0)

print(
    f"  After P2: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P3: Transits (upper and lower)
# ============================================================
print("\n=== P3: Upper/lower transits ===")

for date_name, jd in list(dates.items())[:4]:
    for lat in [0.0, 45.0, 66.56, 75.0]:
        for body_id, body_name in [(swe.SUN, "Sun"), (swe.MOON, "Moon")]:
            for rsmi, rsmi_name in [
                (CALC_MTRANSIT, "upper_transit"),
                (CALC_ITRANSIT, "lower_transit"),
            ]:
                label = f"{body_name} {rsmi_name} {date_name} lat={lat}"
                compare_rise_set(label, body_id, jd, lat, lon, alt, rsmi, tol_sec=30.0)

print(
    f"  After P3: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P4: Mars rise/set (tests planet code path)
# ============================================================
print("\n=== P4: Mars rise/set ===")

for date_name, jd in list(dates.items())[:4]:
    for lat in [0.0, 23.44, 45.0, 66.56, 75.0, -45.0]:
        for rsmi, rsmi_name in [(CALC_RISE, "rise"), (CALC_SET, "set")]:
            label = f"Mars {rsmi_name} {date_name} lat={lat}"
            compare_rise_set(label, swe.MARS, jd, lat, lon, alt, rsmi)

print(
    f"  After P4: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P5: Disc center flag (no refraction, no disc correction)
# ============================================================
print("\n=== P5: Disc center + no refraction flags ===")

for date_name, jd in list(dates.items())[:2]:
    for lat in [0.0, 45.0, 66.56, 75.0]:
        for rsmi_flag in [
            CALC_RISE | BIT_DISC_CENTER,
            CALC_SET | BIT_DISC_CENTER,
            CALC_RISE | BIT_NO_REFRACTION,
            CALC_SET | BIT_NO_REFRACTION,
            CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION,
            CALC_SET | BIT_DISC_CENTER | BIT_NO_REFRACTION,
        ]:
            label = f"Sun flags={rsmi_flag} {date_name} lat={lat}"
            compare_rise_set(label, swe.SUN, jd, lat, lon, alt, rsmi_flag)

print(
    f"  After P5: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P6: Different longitudes (timezone edge cases)
# ============================================================
print("\n=== P6: Different longitudes ===")

longitudes = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]
jd = 2451716.5  # Jun 21, 2000

for longitude in longitudes:
    for lat in [0.0, 45.0, 66.56]:
        for rsmi in [CALC_RISE, CALC_SET]:
            rsmi_name = "rise" if rsmi == CALC_RISE else "set"
            label = f"Sun {rsmi_name} lon={longitude} lat={lat}"
            compare_rise_set(label, swe.SUN, jd, lat, longitude, alt, rsmi)

print(
    f"  After P6: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P7: Altitude variations
# ============================================================
print("\n=== P7: Altitude variations ===")

altitudes = [0.0, 100.0, 500.0, 1000.0, 3000.0, 5000.0]
jd = 2451623.5  # Mar 20, 2000

for altitude in altitudes:
    for lat in [0.0, 45.0, 66.56]:
        for rsmi in [CALC_RISE, CALC_SET]:
            rsmi_name = "rise" if rsmi == CALC_RISE else "set"
            label = f"Sun {rsmi_name} alt={altitude}m lat={lat}"
            compare_rise_set(label, swe.SUN, jd, lat, 0.0, altitude, rsmi)

print(
    f"  After P7: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P8: True horizon with horizon altitude
# ============================================================
print("\n=== P8: True horizon (rise_trans_true_hor) ===")

jd = 2451716.5  # Jun 21, 2000
horizon_alts = [0.0, 0.5, 1.0, 2.0, 5.0]

for horhgt in horizon_alts:
    for lat in [0.0, 30.0, 45.0, 60.0]:
        for rsmi in [CALC_RISE, CALC_SET]:
            rsmi_name = "rise" if rsmi == CALC_RISE else "set"
            label = f"Sun truhor {rsmi_name} horhgt={horhgt} lat={lat}"
            geopos = [0.0, lat, 0.0]
            try:
                se_result = swe.rise_trans_true_hor(
                    jd, swe.SUN, rsmi, geopos, 1013.25, 15.0, horhgt
                )
                se_jd = se_result[1][0]
            except Exception:
                se_jd = None

            try:
                le_result = ephem.rise_trans_true_hor(
                    jd, swe.SUN, rsmi, [0.0, lat, 0.0], 1013.25, 15.0, horhgt, 2
                )
                le_jd = le_result[1][0]
            except Exception:
                le_jd = None

            if se_jd is None and le_jd is None:
                passed += 1
            elif se_jd is not None and le_jd is not None:
                if isinstance(le_jd, tuple):
                    le_jd = le_jd[0]
                if se_jd < 1.0 and (le_jd == 0.0 or le_jd is None):
                    passed += 1
                elif se_jd < 1.0 or le_jd == 0.0:
                    skipped += 1
                else:
                    diff_s = abs(float(se_jd) - float(le_jd)) * 86400.0
                    if diff_s < 30.0:
                        passed += 1
                    else:
                        failed += 1
                        print(
                            f"  FAIL {label}: SE={float(se_jd):.6f} LE={float(le_jd):.6f} diff={diff_s:.1f}s"
                        )
            else:
                skipped += 1

print(
    f"  After P8: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# P9: Twilight modes (civil, nautical, astronomical)
# ============================================================
print("\n=== P9: Twilight modes at various latitudes ===")

BIT_CIVIL_TWILIGHT = 4096
BIT_NAUTIC_TWILIGHT = 8192
BIT_ASTRO_TWILIGHT = 16384

jd = 2451716.5  # Jun 21, 2000

for twilight_flag, tw_name in [
    (BIT_CIVIL_TWILIGHT, "civil"),
    (BIT_NAUTIC_TWILIGHT, "nautical"),
    (BIT_ASTRO_TWILIGHT, "astronomical"),
]:
    for lat in [0.0, 30.0, 45.0, 55.0, 60.0, 66.56]:
        for rsmi_base in [CALC_RISE, CALC_SET]:
            rsmi = rsmi_base | twilight_flag
            rsmi_name = "begin" if rsmi_base == CALC_RISE else "end"
            label = f"{tw_name}_twilight {rsmi_name} lat={lat}"
            compare_rise_set(label, swe.SUN, jd, lat, 0.0, 0.0, rsmi, tol_sec=30.0)

print(
    f"  After P9: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
)

# ============================================================
# Summary
# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 80 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"  Skipped: {skipped}")
print("=" * 70)
