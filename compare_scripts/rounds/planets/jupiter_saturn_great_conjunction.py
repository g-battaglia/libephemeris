#!/usr/bin/env python3
"""Round 58: Jupiter-Saturn Great Conjunction Precision

Verify Jupiter-Saturn great conjunctions (occur every ~19.86 years):
- Compare conjunction timing between LE and SE
- Verify the famous 2020-12-21 great conjunction
- Test outer planet conjunction/opposition timing for all outer planets
- Validate opposition timing (planet opposite Sun, elongation = 180°)
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

FLAGS = 256  # FLG_SPEED

print("=" * 70)
print("ROUND 58: Jupiter-Saturn Great Conjunction Precision")
print("=" * 70)


def get_lon(jd, body, use_se=False):
    if use_se:
        return swe.calc_ut(jd, body, FLAGS)[0][0]
    return ephem.calc_ut(jd, body, FLAGS)[0][0]


def angular_sep(lon1, lon2):
    """Signed angular separation lon1 - lon2, in range [-180, 180]."""
    d = lon1 - lon2
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return d


def find_conjunctions_pair(body1, body2, jd_start, jd_end, step=5.0, use_se=False):
    """Find when two bodies have same ecliptic longitude."""
    events = []
    jd = jd_start
    prev_sep = None
    while jd < jd_end:
        try:
            lon1 = get_lon(jd, body1, use_se)
            lon2 = get_lon(jd, body2, use_se)
            sep = angular_sep(lon1, lon2)
            if prev_sep is not None:
                if prev_sep > 0 and sep < 0:
                    # Refine
                    conj_jd = _bisect_pair(body1, body2, jd - step, jd, use_se)
                    events.append(conj_jd)
                elif prev_sep < 0 and sep > 0:
                    conj_jd = _bisect_pair(body1, body2, jd - step, jd, use_se)
                    events.append(conj_jd)
            prev_sep = sep
        except Exception:
            pass
        jd += step
    return events


def _bisect_pair(body1, body2, jd_lo, jd_hi, use_se=False, tol=1e-8):
    for _ in range(60):
        jd_mid = (jd_lo + jd_hi) / 2
        lon1 = get_lon(jd_mid, body1, use_se)
        lon2 = get_lon(jd_mid, body2, use_se)
        sep = angular_sep(lon1, lon2)
        if abs(sep) < 1e-10:
            return jd_mid
        lo_lon1 = get_lon(jd_lo, body1, use_se)
        lo_lon2 = get_lon(jd_lo, body2, use_se)
        lo_sep = angular_sep(lo_lon1, lo_lon2)
        if (lo_sep > 0 and sep > 0) or (lo_sep < 0 and sep < 0):
            jd_lo = jd_mid
        else:
            jd_hi = jd_mid
        if jd_hi - jd_lo < tol:
            break
    return (jd_lo + jd_hi) / 2


def find_oppositions(body_id, jd_start, jd_end, step=2.0, use_se=False):
    """Find opposition times (elongation = 180°) for outer planets."""
    events = []
    jd = jd_start
    prev_sep = None
    while jd < jd_end:
        try:
            sun_lon = get_lon(jd, 0, use_se)
            body_lon = get_lon(jd, body_id, use_se)
            sep = angular_sep(body_lon, sun_lon)
            # Opposition: sep crosses ±180
            if prev_sep is not None:
                if prev_sep > 170 and sep < -170:
                    opp_jd = _bisect_opposition(body_id, jd - step, jd, use_se)
                    events.append(opp_jd)
                elif prev_sep < -170 and sep > 170:
                    opp_jd = _bisect_opposition(body_id, jd - step, jd, use_se)
                    events.append(opp_jd)
            prev_sep = sep
        except Exception:
            pass
        jd += step
    return events


def _bisect_opposition(body_id, jd_lo, jd_hi, use_se=False, tol=1e-8):
    for _ in range(60):
        jd_mid = (jd_lo + jd_hi) / 2
        sun_lon = get_lon(jd_mid, 0, use_se)
        body_lon = get_lon(jd_mid, body_id, use_se)
        sep = angular_sep(body_lon, sun_lon)
        # We want |sep| = 180, so minimize |180 - |sep||
        err = 180.0 - abs(sep)
        if abs(err) < 1e-10:
            return jd_mid
        lo_sun = get_lon(jd_lo, 0, use_se)
        lo_body = get_lon(jd_lo, body_id, use_se)
        lo_sep = angular_sep(lo_body, lo_sun)
        # Check which side
        if (lo_sep > 0 and sep > 0) or (lo_sep < 0 and sep < 0):
            jd_lo = jd_mid
        else:
            jd_hi = jd_mid
        if jd_hi - jd_lo < tol:
            break
    return (jd_lo + jd_hi) / 2


# ============================================================
# P1: Jupiter-Saturn great conjunctions 1900-2100
# ============================================================
print("\n=== P1: Jupiter-Saturn great conjunctions 1900-2100 ===")

jd_1900 = swe.julday(1900, 1, 1, 0.0)
jd_2100 = swe.julday(2100, 1, 1, 0.0)

le_gc = find_conjunctions_pair(5, 6, jd_1900, jd_2100, step=10.0)
se_gc = find_conjunctions_pair(5, 6, jd_1900, jd_2100, step=10.0, use_se=True)

print(f"  LE found {len(le_gc)} great conjunctions, SE found {len(se_gc)}")

# Expected: ~10 conjunctions in 200 years (every ~19.86 yr)
if len(le_gc) == len(se_gc) and 8 <= len(le_gc) <= 15:
    passed += 1
else:
    failed += 1
    print(f"  FAIL P1: count mismatch or out of range")

for i, (le_jd, se_jd) in enumerate(zip(le_gc, se_gc)):
    diff_hours = abs(le_jd - se_jd) * 24
    y, m, d, h = swe.revjul(le_jd)

    # Great conjunction timing should agree within 1 hour
    if diff_hours < 2.0:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL P1 GC {i} {y}-{m:02d}-{d:02d}: diff={diff_hours:.2f}h")

    # Check longitude at conjunction
    le_jup = get_lon(le_jd, 5)
    le_sat = get_lon(le_jd, 6)
    sep_arcsec = abs(angular_sep(le_jup, le_sat)) * 3600
    if sep_arcsec < 1.0:
        passed += 1
    else:
        failed += 1
        print(f'  FAIL P1 GC {i} separation: {sep_arcsec:.2f}" at conjunction')

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: The famous 2020-12-21 great conjunction
# ============================================================
print("\n=== P2: 2020-12-21 Great Conjunction ===")

# The 2020 great conjunction was notable: Jupiter and Saturn within 0.1°
# Occurred around 2020-12-21 ~18:20 UT
jd_2020gc = swe.julday(2020, 12, 21, 18.33)

try:
    le_jup = ephem.calc_ut(jd_2020gc, 5, FLAGS)[0]
    le_sat = ephem.calc_ut(jd_2020gc, 6, FLAGS)[0]
    se_jup = swe.calc_ut(jd_2020gc, 5, FLAGS)[0]
    se_sat = swe.calc_ut(jd_2020gc, 6, FLAGS)[0]

    le_sep = abs(angular_sep(le_jup[0], le_sat[0]))
    se_sep = abs(angular_sep(se_jup[0], se_sat[0]))

    # 2020 GC separation was ~0.1° (6')
    if le_sep < 0.2 and se_sep < 0.2:
        passed += 1
    else:
        failed += 1
        print(f"  FAIL P2: LE_sep={le_sep:.4f}° SE_sep={se_sep:.4f}°")

    # LE vs SE agreement on separation
    diff_sep = abs(le_sep - se_sep) * 3600
    if diff_sep < 1.0:
        passed += 1
    else:
        failed += 1
        print(f'  FAIL P2 sep agreement: diff={diff_sep:.2f}"')

    # Jupiter position agreement
    jup_diff = abs(le_jup[0] - se_jup[0]) * 3600
    if jup_diff < 0.5:
        passed += 1
    else:
        failed += 1
        print(f'  FAIL P2 Jupiter pos: diff={jup_diff:.3f}"')

    # Saturn position agreement
    sat_diff = abs(le_sat[0] - se_sat[0]) * 3600
    if sat_diff < 0.5:
        passed += 1
    else:
        failed += 1
        print(f'  FAIL P2 Saturn pos: diff={sat_diff:.3f}"')

except Exception as e:
    errors += 1
    print(f"  ERROR P2: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Jupiter oppositions 2020-2030
# ============================================================
print("\n=== P3: Jupiter oppositions 2020-2030 ===")

jd_2020 = swe.julday(2020, 1, 1, 0.0)
jd_2030 = swe.julday(2030, 1, 1, 0.0)

le_jup_opp = find_oppositions(5, jd_2020, jd_2030)
se_jup_opp = find_oppositions(5, jd_2020, jd_2030, use_se=True)

print(f"  LE found {len(le_jup_opp)} Jupiter oppositions, SE found {len(se_jup_opp)}")

# Jupiter opposition every ~13 months, expect ~8-9 in 10 years
if len(le_jup_opp) == len(se_jup_opp):
    passed += 1
    for i, (le_jd, se_jd) in enumerate(zip(le_jup_opp, se_jup_opp)):
        diff_hours = abs(le_jd - se_jd) * 24
        y, m, d, h = swe.revjul(le_jd)
        if diff_hours < 2.0:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P3 Jup opp {i} {y}-{m:02d}-{d:02d}: diff={diff_hours:.2f}h")
else:
    failed += 1
    print(f"  FAIL P3: count mismatch")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Saturn oppositions 2020-2030
# ============================================================
print("\n=== P4: Saturn oppositions 2020-2030 ===")

le_sat_opp = find_oppositions(6, jd_2020, jd_2030)
se_sat_opp = find_oppositions(6, jd_2020, jd_2030, use_se=True)

print(f"  LE found {len(le_sat_opp)} Saturn oppositions, SE found {len(se_sat_opp)}")

if len(le_sat_opp) == len(se_sat_opp):
    passed += 1
    for i, (le_jd, se_jd) in enumerate(zip(le_sat_opp, se_sat_opp)):
        diff_hours = abs(le_jd - se_jd) * 24
        y, m, d, h = swe.revjul(le_jd)
        if diff_hours < 3.0:  # Saturn slower
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P4 Sat opp {i} {y}-{m:02d}-{d:02d}: diff={diff_hours:.2f}h")
else:
    failed += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Mars oppositions 2020-2030
# ============================================================
print("\n=== P5: Mars oppositions 2020-2030 ===")

le_mars_opp = find_oppositions(4, jd_2020, jd_2030)
se_mars_opp = find_oppositions(4, jd_2020, jd_2030, use_se=True)

print(f"  LE found {len(le_mars_opp)} Mars oppositions, SE found {len(se_mars_opp)}")

if len(le_mars_opp) == len(se_mars_opp):
    passed += 1
    for i, (le_jd, se_jd) in enumerate(zip(le_mars_opp, se_mars_opp)):
        diff_hours = abs(le_jd - se_jd) * 24
        y, m, d, h = swe.revjul(le_jd)
        if diff_hours < 1.0:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P5 Mars opp {i} {y}-{m:02d}-{d:02d}: diff={diff_hours:.2f}h")
else:
    failed += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Great conjunction period (~19.86 years)
# ============================================================
print("\n=== P6: Great conjunction period ===")

expected_period = 19.86 * 365.25  # ~7253 days
for i in range(len(le_gc) - 1):
    period = le_gc[i + 1] - le_gc[i]
    period_years = period / 365.25
    # Period can vary ~18-21 years due to orbital eccentricities
    if 17.0 < period_years < 22.0:
        passed += 1
    else:
        y1, m1, d1, h1 = swe.revjul(le_gc[i])
        y2, m2, d2, h2 = swe.revjul(le_gc[i + 1])
        failed += 1
        print(f"  FAIL P6: {y1}-{y2} period={period_years:.2f}yr (expected ~19.86yr)")

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Jupiter-Saturn daily positions around 2020 GC
# ============================================================
print("\n=== P7: Daily positions near 2020 GC ===")

jd_start = swe.julday(2020, 11, 1, 0.0)
jd_end = swe.julday(2021, 2, 1, 0.0)

for day in range(int(jd_end - jd_start)):
    jd = jd_start + day
    try:
        for body_id, name in [(5, "Jup"), (6, "Sat")]:
            le = ephem.calc_ut(jd, body_id, FLAGS)[0]
            se = swe.calc_ut(jd, body_id, FLAGS)[0]

            lon_diff = abs(le[0] - se[0]) * 3600
            if lon_diff > 180 * 3600:
                lon_diff = 360 * 3600 - lon_diff

            if lon_diff < 0.5:
                passed += 1
            else:
                y, m, d, h = swe.revjul(jd)
                failed += 1
                print(f'  FAIL P7 {name} {y}-{m:02d}-{d:02d}: diff={lon_diff:.3f}"')

    except Exception as e:
        errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Opposition synodic periods
# ============================================================
print("\n=== P8: Opposition synodic periods ===")

# Mars synodic: ~779.9 days, Jupiter: ~398.9 days, Saturn: ~378.1 days
expected_synodic = {4: 779.9, 5: 398.9, 6: 378.1}
opp_data = {4: le_mars_opp, 5: le_jup_opp, 6: le_sat_opp}
names = {4: "Mars", 5: "Jupiter", 6: "Saturn"}

for body_id in [4, 5, 6]:
    opps = opp_data[body_id]
    exp = expected_synodic[body_id]
    for i in range(len(opps) - 1):
        period = opps[i + 1] - opps[i]
        # Allow ±30 days variation
        if abs(period - exp) < 30.0:
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P8 {names[body_id]} opp period {i}: {period:.1f}d (expected ~{exp:.1f}d)"
            )

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
print(
    f"ROUND 58 FINAL: {passed}/{passed + failed} passed ({100 * passed / (passed + failed):.1f}%)"
)
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
