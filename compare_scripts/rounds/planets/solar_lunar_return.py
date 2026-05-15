#!/usr/bin/env python3
"""Round 89: Solar and Lunar Return Timing

Tests precise timing of solar returns (Sun returning to natal position)
and lunar returns (Moon returning to natal position) using solcross_ut
and mooncross_ut. Verifies sub-second precision.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0

print("=" * 70)
print("ROUND 89: Solar and Lunar Return Timing")
print("=" * 70)

# ============================================================
# P1: Solar returns — Sun crossing natal longitude each year
# ============================================================
print("\n=== P1: Solar returns (10 natal charts x 20 years) ===")

# 10 natal Sun positions (various dates)
natal_dates = [
    (1970, 3, 21),
    (1975, 7, 4),
    (1980, 1, 15),
    (1985, 10, 31),
    (1990, 6, 21),
    (1960, 12, 25),
    (1955, 4, 10),
    (1965, 8, 15),
    (1972, 11, 11),
    (1988, 2, 29),
]

for ny, nm, nd in natal_dates:
    natal_jd = swe.julday(ny, nm, nd, 12.0)
    se_natal = swe.calc_ut(natal_jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED)
    natal_lon = se_natal[0][0]

    for year_offset in range(1, 21):
        search_jd = natal_jd + year_offset * 365.0
        label = f"solar_return {ny}-{nm:02d}-{nd:02d} +{year_offset}yr"
        try:
            se_jd = swe.solcross_ut(natal_lon, search_jd, swe.FLG_SWIEPH)
            le_jd = ephem.solcross_ut(natal_lon, search_jd, 2)
            diff_s = abs(se_jd - le_jd) * 86400.0
            if diff_s < 2.0:  # 2 second tolerance
                passed += 1
            else:
                failed += 1
                print(f"  FAIL {label}: diff={diff_s:.2f}s")
        except Exception as e:
            errors += 1

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Lunar returns — Moon crossing natal longitude each month
# ============================================================
print("\n=== P2: Lunar returns (5 natal charts x 24 months) ===")

for ny, nm, nd in natal_dates[:5]:
    natal_jd = swe.julday(ny, nm, nd, 12.0)
    se_natal = swe.calc_ut(natal_jd, swe.MOON, swe.FLG_SWIEPH | swe.FLG_SPEED)
    natal_moon_lon = se_natal[0][0]

    for month_offset in range(1, 25):
        search_jd = natal_jd + month_offset * 27.3  # ~sidereal month
        label = f"lunar_return {ny}-{nm:02d}-{nd:02d} +{month_offset}m"
        try:
            se_jd = swe.mooncross_ut(natal_moon_lon, search_jd, swe.FLG_SWIEPH)
            le_jd = ephem.mooncross_ut(natal_moon_lon, search_jd, 2)
            diff_s = abs(se_jd - le_jd) * 86400.0
            if diff_s < 1.0:  # 1 second tolerance (Moon moves fast)
                passed += 1
            else:
                failed += 1
                print(f"  FAIL {label}: diff={diff_s:.2f}s")
        except Exception as e:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Solar return position verification
# ============================================================
print("\n=== P3: Solar return position at crossing ===")

for ny, nm, nd in natal_dates[:5]:
    natal_jd = swe.julday(ny, nm, nd, 12.0)
    se_natal = swe.calc_ut(natal_jd, swe.SUN, swe.FLG_SWIEPH)
    natal_lon = se_natal[0][0]

    for year_offset in [1, 5, 10, 20]:
        search_jd = natal_jd + year_offset * 365.0
        label = f"sr_pos {ny} +{year_offset}yr"
        try:
            le_jd = ephem.solcross_ut(natal_lon, search_jd, 2)
            le_pos = ephem.calc_ut(le_jd, 0, 2 | 256)
            actual_lon = le_pos[0][0]
            diff_arcsec = abs(actual_lon - natal_lon) * 3600.0
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec
            if diff_arcsec < 0.001:  # sub-milliarcsecond
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: target={natal_lon:.6f} actual={actual_lon:.6f} diff={diff_arcsec:.6f}"'
                )
        except Exception as e:
            errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Lunar return position verification
# ============================================================
print("\n=== P4: Lunar return position at crossing ===")

for ny, nm, nd in natal_dates[:5]:
    natal_jd = swe.julday(ny, nm, nd, 12.0)
    se_natal = swe.calc_ut(natal_jd, swe.MOON, swe.FLG_SWIEPH)
    natal_moon_lon = se_natal[0][0]

    for month_offset in [1, 6, 12]:
        search_jd = natal_jd + month_offset * 27.3
        label = f"lr_pos {ny} +{month_offset}m"
        try:
            le_jd = ephem.mooncross_ut(natal_moon_lon, search_jd, 2)
            le_pos = ephem.calc_ut(le_jd, 1, 2 | 256)
            actual_lon = le_pos[0][0]
            diff_arcsec = abs(actual_lon - natal_moon_lon) * 3600.0
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec
            if diff_arcsec < 0.01:  # 10 mas
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: target={natal_moon_lon:.6f} actual={actual_lon:.6f} diff={diff_arcsec:.6f}"'
                )
        except Exception as e:
            errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Consecutive solar returns — verify ~365.25 day spacing
# ============================================================
print("\n=== P5: Solar return interval consistency ===")

natal_jd = swe.julday(1980, 1, 1, 12.0)
se_natal = swe.calc_ut(natal_jd, swe.SUN, swe.FLG_SWIEPH)
natal_lon = se_natal[0][0]

prev_jd = None
for year_offset in range(1, 31):
    search_jd = natal_jd + year_offset * 365.0
    label = f"sr_interval +{year_offset}yr"
    try:
        le_jd = ephem.solcross_ut(natal_lon, search_jd, 2)
        if prev_jd is not None:
            interval = le_jd - prev_jd
            # Tropical year is ~365.2422 days, but varies slightly
            if 365.0 < interval < 366.0:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL {label}: interval={interval:.4f} days")
        prev_jd = le_jd
    except Exception as e:
        errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 89 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
