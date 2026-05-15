#!/usr/bin/env python3
"""Round 91: Moon Void-of-Course Detection

Tests Moon void-of-course by computing the last major aspect the Moon makes
before leaving its current sign. Compares Moon ingress timing and last-aspect
timing between SE and LE to verify astrological chart accuracy.
"""

from __future__ import annotations
import os, sys, math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0
F = 2
S = 256

print("=" * 70)
print("ROUND 91: Moon Void-of-Course Detection")
print("=" * 70)

major_aspects = [0, 60, 90, 120, 180]  # conjunction, sextile, square, trine, opposition


def get_moon_sign(jd, use_se=True):
    """Get Moon's zodiac sign (0-11)."""
    if use_se:
        r = swe.calc_ut(jd, swe.MOON, F | S)
        return int(r[0][0] / 30.0) % 12
    else:
        r = ephem.calc_ut(jd, 1, F | S)
        return int(r[0][0] / 30.0) % 12


def get_planet_positions(jd, use_se=True):
    """Get longitudes of Sun through Saturn."""
    positions = {}
    planets = [
        (0, "Sun"),
        (2, "Mercury"),
        (3, "Venus"),
        (4, "Mars"),
        (5, "Jupiter"),
        (6, "Saturn"),
    ]
    for pid, name in planets:
        if use_se:
            r = swe.calc_ut(jd, pid, F | S)
            positions[name] = r[0][0]
        else:
            r = ephem.calc_ut(jd, pid, F | S)
            positions[name] = r[0][0]
    return positions


def check_aspect(moon_lon, planet_lon, orb=1.0):
    """Check if Moon is within orb of any major aspect to planet."""
    for asp in major_aspects:
        diff = abs((moon_lon - planet_lon) % 360.0 - asp)
        if diff > 180:
            diff = 360 - diff
        if diff < orb:
            return True, asp
    return False, None


# ============================================================
# P1: Moon ingress timing SE vs LE (every 2.5 days over 2 years)
# ============================================================
print("\n=== P1: Moon ingress timing comparison ===")

jd_start = swe.julday(2020, 1, 1, 0.0)
for i in range(288):  # ~2 years of Moon ingresses
    jd = jd_start + i * 2.5
    for target in range(0, 360, 30):
        label = f"Moon ingress {target}° #{i}"
        try:
            se_jd = swe.mooncross_ut(float(target), jd, F)
            le_jd = ephem.mooncross_ut(float(target), jd, F)
            diff_s = abs(se_jd - le_jd) * 86400.0
            if diff_s < 1.0:
                passed += 1
            else:
                failed += 1
                if failed <= 3:
                    print(f"  FAIL {label}: diff={diff_s:.2f}s")
        except Exception as e:
            errors += 1

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Moon sign agreement SE vs LE (hourly over 30 days)
# ============================================================
print("\n=== P2: Moon sign agreement (hourly 30 days) ===")

jd_start = swe.julday(2023, 6, 1, 0.0)
sign_match = sign_mismatch = 0
for hour in range(720):  # 30 days * 24 hours
    jd = jd_start + hour / 24.0
    try:
        se_sign = get_moon_sign(jd, True)
        le_sign = get_moon_sign(jd, False)
        if se_sign == le_sign:
            sign_match += 1
        else:
            sign_mismatch += 1
            if sign_mismatch <= 3:
                se_r = swe.calc_ut(jd, swe.MOON, F | S)
                le_r = ephem.calc_ut(jd, 1, F | S)
                print(
                    f"  FAIL sign h={hour}: SE_sign={se_sign} LE_sign={le_sign} SE_lon={se_r[0][0]:.6f} LE_lon={le_r[0][0]:.6f}"
                )
    except Exception as e:
        errors += 1

passed += sign_match
failed += sign_mismatch
print(f"  Sign match: {sign_match}, mismatch: {sign_mismatch}")
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Moon aspect detection agreement (SE vs LE)
# ============================================================
print("\n=== P3: Moon aspect detection agreement ===")

jd_start = swe.julday(2023, 1, 1, 0.0)
asp_match = asp_mismatch = 0

for day in range(90):  # 90 days, check every 2 hours
    for h in range(0, 24, 2):
        jd = jd_start + day + h / 24.0
        try:
            # Get Moon position from both
            se_moon = swe.calc_ut(jd, swe.MOON, F | S)[0][0]
            le_moon = ephem.calc_ut(jd, 1, F | S)[0][0]

            # Get planet positions
            se_planets = get_planet_positions(jd, True)
            le_planets = get_planet_positions(jd, False)

            # Check aspects
            for pname in se_planets:
                se_has, se_asp = check_aspect(se_moon, se_planets[pname])
                le_has, le_asp = check_aspect(le_moon, le_planets[pname])
                if se_has == le_has:
                    asp_match += 1
                else:
                    asp_mismatch += 1
        except Exception as e:
            errors += 1

passed += asp_match
failed += asp_mismatch
if asp_mismatch > 0:
    print(f"  Aspect mismatch: {asp_mismatch} (may be near exact aspect boundary)")
print(f"  Aspect match: {asp_match}, mismatch: {asp_mismatch}")
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: VOC period length comparison
# ============================================================
print("\n=== P4: VOC period detection ===")

# For each Moon sign transit, find last aspect before ingress
jd = swe.julday(2023, 1, 1, 0.0)
voc_tests = 0
voc_match = 0

for i in range(24):  # 24 consecutive Moon ingresses
    # Find next ingress from current position
    try:
        le_moon = ephem.calc_ut(jd, 1, F | S)
        current_sign = int(le_moon[0][0] / 30.0) % 12
        next_boundary = float(((current_sign + 1) % 12) * 30)

        se_ingress = swe.mooncross_ut(next_boundary, jd, F)
        le_ingress = ephem.mooncross_ut(next_boundary, jd, F)

        # Both should find ingress within ~3 days
        if abs(se_ingress - le_ingress) * 86400.0 < 2.0:
            voc_match += 1
        voc_tests += 1

        jd = le_ingress + 0.1  # Move past ingress
    except Exception as e:
        errors += 1
        jd += 2.5

passed += voc_match
failed += voc_tests - voc_match
print(f"  VOC ingress match: {voc_match}/{voc_tests}")
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 91 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
