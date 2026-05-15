#!/usr/bin/env python3
"""Round 57: Venus Inferior/Superior Conjunction Timing

Verify Venus conjunction events:
- Inferior conjunction: Venus between Earth and Sun (elongation ~0°, Venus retrograde)
- Superior conjunction: Venus on far side of Sun (elongation ~0°, Venus direct)
- Venus synodic period: ~583.9 days
- Also test elongation precision for all inner planets
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

SUN = 0
MERCURY = 2
VENUS = 3
FLAGS = 256  # FLG_SPEED

print("=" * 70)
print("ROUND 57: Venus Inferior/Superior Conjunction Timing")
print("=" * 70)


def elongation(jd, body_id, use_se=False):
    """Calculate elongation of a body from the Sun."""
    if use_se:
        sun = swe.calc_ut(jd, SUN, FLAGS)[0]
        body = swe.calc_ut(jd, body_id, FLAGS)[0]
    else:
        sun = ephem.calc_ut(jd, SUN, FLAGS)[0]
        body = ephem.calc_ut(jd, body_id, FLAGS)[0]

    diff = body[0] - sun[0]
    if diff > 180:
        diff -= 360
    if diff < -180:
        diff += 360
    return diff


def find_conjunctions(body_id, jd_start, jd_end, step=1.0):
    """Find conjunction times (elongation crosses 0)."""
    conjunctions = []
    jd = jd_start
    prev_elong = None
    while jd < jd_end:
        try:
            elong = elongation(jd, body_id)
            if prev_elong is not None:
                # Sign change = conjunction
                if prev_elong > 0 and elong < 0:
                    # Crossed 0 going negative
                    conj_jd = _bisect_conjunction(body_id, jd - step, jd)
                    # Check if inferior (retrograde) or superior (direct)
                    result = ephem.calc_ut(conj_jd, body_id, FLAGS)
                    speed = result[0][3]
                    ctype = "inf" if speed < 0 else "sup"
                    conjunctions.append((ctype, conj_jd))
                elif prev_elong < 0 and elong > 0:
                    conj_jd = _bisect_conjunction(body_id, jd - step, jd)
                    result = ephem.calc_ut(conj_jd, body_id, FLAGS)
                    speed = result[0][3]
                    ctype = "inf" if speed < 0 else "sup"
                    conjunctions.append((ctype, conj_jd))
            prev_elong = elong
        except Exception:
            pass
        jd += step
    return conjunctions


def _bisect_conjunction(body_id, jd_lo, jd_hi, tol=1e-8):
    """Bisect to find when elongation = 0."""
    for _ in range(60):
        jd_mid = (jd_lo + jd_hi) / 2
        elong = elongation(jd_mid, body_id)
        if abs(elong) < 1e-10:
            return jd_mid
        lo_elong = elongation(jd_lo, body_id)
        if (lo_elong > 0 and elong > 0) or (lo_elong < 0 and elong < 0):
            jd_lo = jd_mid
        else:
            jd_hi = jd_mid
        if jd_hi - jd_lo < tol:
            break
    return (jd_lo + jd_hi) / 2


def find_conjunctions_se(body_id, jd_start, jd_end, step=1.0):
    """Find conjunction times using SE."""
    conjunctions = []
    jd = jd_start
    prev_elong = None
    while jd < jd_end:
        try:
            elong = elongation(jd, body_id, use_se=True)
            if prev_elong is not None:
                if (prev_elong > 0 and elong < 0) or (prev_elong < 0 and elong > 0):
                    conj_jd = _bisect_conjunction_se(body_id, jd - step, jd)
                    result = swe.calc_ut(conj_jd, body_id, FLAGS)
                    speed = result[0][3]
                    ctype = "inf" if speed < 0 else "sup"
                    conjunctions.append((ctype, conj_jd))
            prev_elong = elong
        except Exception:
            pass
        jd += step
    return conjunctions


def _bisect_conjunction_se(body_id, jd_lo, jd_hi, tol=1e-8):
    for _ in range(60):
        jd_mid = (jd_lo + jd_hi) / 2
        elong = elongation(jd_mid, body_id, use_se=True)
        if abs(elong) < 1e-10:
            return jd_mid
        lo_elong = elongation(jd_lo, body_id, use_se=True)
        if (lo_elong > 0 and elong > 0) or (lo_elong < 0 and elong < 0):
            jd_lo = jd_mid
        else:
            jd_hi = jd_mid
        if jd_hi - jd_lo < tol:
            break
    return (jd_lo + jd_hi) / 2


# ============================================================
# P1: Venus conjunctions 2020-2030
# ============================================================
print("\n=== P1: Venus conjunction timing 2020-2030 ===")

jd_2020 = swe.julday(2020, 1, 1, 0.0)
jd_2030 = swe.julday(2030, 1, 1, 0.0)

le_venus = find_conjunctions(VENUS, jd_2020, jd_2030, step=0.5)
se_venus = find_conjunctions_se(VENUS, jd_2020, jd_2030, step=0.5)

print(f"  LE found {len(le_venus)} Venus conjunctions, SE found {len(se_venus)}")

if len(le_venus) == len(se_venus):
    passed += 1
    for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(zip(le_venus, se_venus)):
        if le_type == se_type:
            diff_sec = abs(le_jd - se_jd) * 86400
            if diff_sec < 60.0:
                passed += 1
            else:
                failed += 1
                y, m, d, h = swe.revjul(le_jd)
                print(
                    f"  FAIL P1 conj {i} ({le_type}) {y}-{m:02d}-{d:02d}: diff={diff_sec:.1f}s"
                )
        else:
            failed += 1
            print(f"  FAIL P1 conj {i}: type mismatch LE={le_type} SE={se_type}")
else:
    failed += 1
    print(f"  FAIL P1: count mismatch")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P2: Venus synodic period verification
# ============================================================
print("\n=== P2: Venus synodic period ===")

# Expected synodic period: 583.9 days
expected_synodic = 583.9

for i in range(len(le_venus) - 1):
    le_type1, le_jd1 = le_venus[i]
    le_type2, le_jd2 = le_venus[i + 1]
    if le_type1 == le_type2:
        # Same type conjunction = one synodic period
        period = le_jd2 - le_jd1
        if abs(period - expected_synodic) < 5.0:  # Within 5 days of expected
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P2 period {i}: {period:.1f}d (expected ~{expected_synodic:.1f}d)"
            )

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P3: Mercury conjunctions 2020-2030
# ============================================================
print("\n=== P3: Mercury conjunction timing 2020-2030 ===")

le_merc = find_conjunctions(MERCURY, jd_2020, jd_2030, step=0.5)
se_merc = find_conjunctions_se(MERCURY, jd_2020, jd_2030, step=0.5)

print(f"  LE found {len(le_merc)} Mercury conjunctions, SE found {len(se_merc)}")

if len(le_merc) == len(se_merc):
    passed += 1
    for i, ((le_type, le_jd), (se_type, se_jd)) in enumerate(zip(le_merc, se_merc)):
        if le_type == se_type:
            diff_sec = abs(le_jd - se_jd) * 86400
            if diff_sec < 30.0:
                passed += 1
            else:
                failed += 1
                y, m, d, h = swe.revjul(le_jd)
                print(
                    f"  FAIL P3 conj {i} ({le_type}) {y}-{m:02d}-{d:02d}: diff={diff_sec:.1f}s"
                )
        else:
            failed += 1
else:
    failed += 1
    print(f"  FAIL P3: count mismatch LE={len(le_merc)} SE={len(se_merc)}")

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P4: Mercury synodic period
# ============================================================
print("\n=== P4: Mercury synodic period ===")

# Expected Mercury synodic period: 115.9 days
expected_merc_synodic = 115.9

for i in range(len(le_merc) - 1):
    le_type1, le_jd1 = le_merc[i]
    le_type2, le_jd2 = le_merc[i + 1]
    if le_type1 == le_type2:
        period = le_jd2 - le_jd1
        if abs(period - expected_merc_synodic) < 15.0:  # Mercury has more variation
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL P4 period {i}: {period:.1f}d (expected ~{expected_merc_synodic:.1f}d)"
            )

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P5: Elongation precision at various points
# ============================================================
print("\n=== P5: Venus elongation precision (daily for 2024) ===")

jd_start = swe.julday(2024, 1, 1, 0.0)
for day in range(365):
    jd = jd_start + day
    try:
        le_elong = elongation(jd, VENUS)
        se_elong = elongation(jd, VENUS, use_se=True)

        diff_arcsec = abs(le_elong - se_elong) * 3600
        if diff_arcsec < 0.5:
            passed += 1
        else:
            y, m, d, h = swe.revjul(jd)
            failed += 1
            print(
                f'  FAIL P5 {y}-{m:02d}-{d:02d}: LE={le_elong:.4f}° SE={se_elong:.4f}° diff={diff_arcsec:.2f}"'
            )

    except Exception as e:
        errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P6: Mercury elongation precision (daily for 2024)
# ============================================================
print("\n=== P6: Mercury elongation precision (daily for 2024) ===")

for day in range(365):
    jd = jd_start + day
    try:
        le_elong = elongation(jd, MERCURY)
        se_elong = elongation(jd, MERCURY, use_se=True)

        diff_arcsec = abs(le_elong - se_elong) * 3600
        if diff_arcsec < 0.5:
            passed += 1
        else:
            y, m, d, h = swe.revjul(jd)
            failed += 1
            print(
                f'  FAIL P6 {y}-{m:02d}-{d:02d}: LE={le_elong:.4f}° SE={se_elong:.4f}° diff={diff_arcsec:.2f}"'
            )

    except Exception as e:
        errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P7: Venus maximum elongation values
# ============================================================
print("\n=== P7: Venus max elongation (should be ~46°) ===")

# Find max elongation events for Venus in 2020-2030
max_elongs = []
jd = jd_2020
prev_elong = None
prev_abs_elong = None
while jd < jd_2030:
    try:
        elong = elongation(jd, VENUS)
        abs_elong = abs(elong)
        if prev_abs_elong is not None and prev_elong is not None:
            # Local max when derivative changes sign
            prev_abs2 = abs(prev_elong)
            if prev_abs_elong > abs_elong and prev_abs_elong > prev_abs2:
                max_elongs.append((jd - 1.0, prev_elong))
        prev_abs2 = prev_abs_elong
        prev_abs_elong = abs_elong
        prev_elong = elong
    except Exception:
        pass
    jd += 1.0

for jd_max, elong_val in max_elongs:
    abs_elong = abs(elong_val)
    # Venus max elongation should be 45°-47.5°
    if 44.0 < abs_elong < 48.0:
        passed += 1
    else:
        y, m, d, h = swe.revjul(jd_max)
        failed += 1
        print(
            f"  FAIL P7 {y}-{m:02d}-{d:02d}: max_elong={abs_elong:.2f}° (expected 45-47.5°)"
        )

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
# P8: Venus position at conjunction (should match Sun longitude)
# ============================================================
print("\n=== P8: Venus longitude at conjunction ===")

for i, (ctype, conj_jd) in enumerate(le_venus):
    try:
        le_sun = ephem.calc_ut(conj_jd, SUN, FLAGS)[0]
        le_venus_pos = ephem.calc_ut(conj_jd, VENUS, FLAGS)[0]

        diff = abs(le_venus_pos[0] - le_sun[0])
        if diff > 180:
            diff = 360 - diff

        # At conjunction, Venus should be within ~0.01° of Sun
        if diff < 0.01:
            passed += 1
        else:
            y, m, d, h = swe.revjul(conj_jd)
            failed += 1
            print(
                f"  FAIL P8 conj {i} ({ctype}) {y}-{m:02d}: Venus-Sun diff={diff:.4f}°"
            )

    except Exception as e:
        errors += 1

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")


# ============================================================
print(f"\n{'=' * 70}")
print(
    f"ROUND 57 FINAL: {passed}/{passed + failed} passed ({100 * passed / (passed + failed):.1f}%)"
)
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")

if failed > 0:
    print(f"\n--- FAILURES ({failed}) ---")
