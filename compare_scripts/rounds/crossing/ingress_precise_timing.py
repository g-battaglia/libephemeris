#!/usr/bin/env python3
"""Round 90: Planetary Ingress Precise Timing

Tests precise timing of planetary ingresses (planets crossing sign boundaries
at 0°, 30°, 60°, etc.) using solcross_ut for Sun and cross_ut equivalent
via mooncross_ut for Moon. Compares SE vs LE timing.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0

print("=" * 70)
print("ROUND 90: Planetary Ingress Precise Timing")
print("=" * 70)

sign_boundaries = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
sign_names = [
    "Ari",
    "Tau",
    "Gem",
    "Can",
    "Leo",
    "Vir",
    "Lib",
    "Sco",
    "Sag",
    "Cap",
    "Aqu",
    "Pis",
]

# ============================================================
# P1: Sun ingresses (all 12 signs x 25 years)
# ============================================================
print("\n=== P1: Sun ingresses (12 signs x 25 years) ===")

for year in range(1980, 2026):
    jd_start = swe.julday(year, 1, 1, 0.0)
    for i, target in enumerate(sign_boundaries):
        label = f"Sun->{sign_names[i]} {year}"
        try:
            se_jd = swe.solcross_ut(float(target), jd_start, swe.FLG_SWIEPH)
            le_jd = ephem.solcross_ut(float(target), jd_start, 2)
            diff_s = abs(se_jd - le_jd) * 86400.0
            if diff_s < 2.0:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(f"  FAIL {label}: diff={diff_s:.2f}s")
        except Exception as e:
            errors += 1

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Moon ingresses (12 signs x 24 months)
# ============================================================
print("\n=== P2: Moon ingresses (12 signs x 24 months) ===")

for month in range(24):
    jd_start = swe.julday(2000, 1, 1, 0.0) + month * 29.53
    for i, target in enumerate(sign_boundaries):
        label = f"Moon->{sign_names[i]} m={month}"
        try:
            se_jd = swe.mooncross_ut(float(target), jd_start, swe.FLG_SWIEPH)
            le_jd = ephem.mooncross_ut(float(target), jd_start, 2)
            diff_s = abs(se_jd - le_jd) * 86400.0
            if diff_s < 1.0:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(f"  FAIL {label}: diff={diff_s:.2f}s")
        except Exception as e:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Sun ingress position verification (is Sun really at boundary?)
# ============================================================
print("\n=== P3: Sun at sign boundary verification ===")

for year in range(2000, 2026):
    jd_start = swe.julday(year, 1, 1, 0.0)
    for target in [0, 90, 180, 270]:
        label = f"Sun@{target}° {year}"
        try:
            le_jd = ephem.solcross_ut(float(target), jd_start, 2)
            le_pos = ephem.calc_ut(le_jd, 0, 2 | 256)
            actual = le_pos[0][0]
            diff = abs(actual - target) * 3600.0
            if diff > 180 * 3600:
                diff = 360 * 3600 - diff
            if diff < 0.001:
                passed += 1
            else:
                failed += 1
                print(f'  FAIL {label}: actual={actual:.8f} diff={diff:.6f}"')
        except Exception as e:
            errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Moon ingress position verification
# ============================================================
print("\n=== P4: Moon at sign boundary verification ===")

for month in range(12):
    jd_start = swe.julday(2020, 1, 1, 0.0) + month * 29.53
    for target in [0, 90, 180, 270]:
        label = f"Moon@{target}° m={month}"
        try:
            le_jd = ephem.mooncross_ut(float(target), jd_start, 2)
            le_pos = ephem.calc_ut(le_jd, 1, 2 | 256)
            actual = le_pos[0][0]
            diff = abs(actual - target) * 3600.0
            if diff > 180 * 3600:
                diff = 360 * 3600 - diff
            if diff < 0.01:
                passed += 1
            else:
                failed += 1
                print(f'  FAIL {label}: actual={actual:.8f} diff={diff:.6f}"')
        except Exception as e:
            errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Cardinal ingress timing (equinoxes/solstices) vs SE
# ============================================================
print("\n=== P5: Cardinal ingresses (equinoxes/solstices) ===")

for year in range(1950, 2051):
    jd_start = swe.julday(year, 1, 1, 0.0)
    for target, name in [(0, "VE"), (90, "SS"), (180, "AE"), (270, "WS")]:
        label = f"{name} {year}"
        try:
            se_jd = swe.solcross_ut(float(target), jd_start, swe.FLG_SWIEPH)
            le_jd = ephem.solcross_ut(float(target), jd_start, 2)
            diff_s = abs(se_jd - le_jd) * 86400.0
            if diff_s < 2.0:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(f"  FAIL {label}: diff={diff_s:.2f}s")
        except Exception as e:
            errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: Consecutive Moon ingresses — verify spacing is ~2.3 days
# ============================================================
print("\n=== P6: Moon ingress spacing consistency ===")

jd = swe.julday(2020, 1, 1, 0.0)
prev_jd = None
for i in range(36):  # 36 consecutive ingresses = ~3 months
    target = (i * 30) % 360
    label = f"Moon spacing #{i}"
    try:
        le_jd = ephem.mooncross_ut(float(target), jd, 2)
        if prev_jd is not None:
            spacing = le_jd - prev_jd
            if 1.5 < spacing < 3.5:  # Moon spends ~2.3 days per sign
                passed += 1
            else:
                failed += 1
                print(f"  FAIL {label}: spacing={spacing:.4f} days")
        prev_jd = le_jd
        jd = le_jd + 0.1  # Start searching just after
    except Exception as e:
        errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 90 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
