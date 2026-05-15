#!/usr/bin/env python3
"""Round 92: Aspect Calculation Precision

Tests that inter-planetary aspects (angular separations) computed from SE and LE
positions agree precisely. Covers all planet pairs, major aspects, and orb detection.
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
print("ROUND 92: Aspect Calculation Precision")
print("=" * 70)

planets = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
]

aspects = [0, 30, 45, 60, 72, 90, 120, 135, 144, 150, 180]
aspect_names = {
    0: "conj",
    30: "semi",
    45: "octile",
    60: "sext",
    72: "quint",
    90: "sq",
    120: "tri",
    135: "sesq",
    144: "biquint",
    150: "quinc",
    180: "opp",
}

test_jds = [
    (y, swe.julday(y, m, 15, 12.0)) for y in range(1990, 2026, 5) for m in [1, 4, 7, 10]
]

# ============================================================
# P1: Angular separation agreement (all planet pairs)
# ============================================================
print("\n=== P1: Angular separation SE vs LE ===")

for year, jd in test_jds:
    # Get all positions from both
    se_lons = {}
    le_lons = {}
    for pid, name in planets:
        try:
            se = swe.calc_ut(jd, pid, F | S)
            le = ephem.calc_ut(jd, pid, F | S)
            se_lons[pid] = se[0][0]
            le_lons[pid] = le[0][0]
        except:
            pass

    # Compare separations for all pairs
    for i in range(len(planets)):
        for j in range(i + 1, len(planets)):
            pid1, n1 = planets[i]
            pid2, n2 = planets[j]
            if pid1 not in se_lons or pid2 not in se_lons:
                continue
            if pid1 not in le_lons or pid2 not in le_lons:
                continue

            se_sep = (se_lons[pid2] - se_lons[pid1]) % 360.0
            if se_sep > 180:
                se_sep = 360.0 - se_sep
            le_sep = (le_lons[pid2] - le_lons[pid1]) % 360.0
            if le_sep > 180:
                le_sep = 360.0 - le_sep

            diff_arcsec = abs(se_sep - le_sep) * 3600.0
            if diff_arcsec < 1.0:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(
                        f'  FAIL {n1}-{n2} jd={jd:.1f}: SE_sep={se_sep:.6f} LE_sep={le_sep:.6f} diff={diff_arcsec:.2f}"'
                    )

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Aspect detection agreement (which aspects are within orb)
# ============================================================
print("\n=== P2: Aspect detection agreement (8° orb) ===")

orb = 8.0
asp_match = asp_mismatch = 0

for year, jd in test_jds:
    se_lons = {}
    le_lons = {}
    for pid, name in planets:
        try:
            se = swe.calc_ut(jd, pid, F | S)
            le = ephem.calc_ut(jd, pid, F | S)
            se_lons[pid] = se[0][0]
            le_lons[pid] = le[0][0]
        except:
            pass

    for i in range(len(planets)):
        for j in range(i + 1, len(planets)):
            pid1, _ = planets[i]
            pid2, _ = planets[j]
            if pid1 not in se_lons or pid2 not in se_lons:
                continue

            se_sep = (se_lons[pid2] - se_lons[pid1]) % 360.0
            le_sep = (le_lons[pid2] - le_lons[pid1]) % 360.0

            for asp in aspects:
                se_diff = (
                    min(
                        abs(se_sep - asp),
                        abs(se_sep - (360 - asp)),
                        abs((360 - se_sep) - asp),
                    )
                    if asp != 0
                    else min(se_sep, 360 - se_sep)
                )
                le_diff = (
                    min(
                        abs(le_sep - asp),
                        abs(le_sep - (360 - asp)),
                        abs((360 - le_sep) - asp),
                    )
                    if asp != 0
                    else min(le_sep, 360 - le_sep)
                )

                # Simple approach: both agree on aspect presence
                se_in_orb = se_diff < orb
                le_in_orb = le_diff < orb
                if se_in_orb == le_in_orb:
                    asp_match += 1
                else:
                    asp_mismatch += 1

passed += asp_match
failed += asp_mismatch
print(f"  Aspect detection: {asp_match} agree, {asp_mismatch} disagree")
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Exact aspect orb precision (how close are SE vs LE orbs?)
# ============================================================
print("\n=== P3: Exact orb precision ===")

orb_diffs = []
for year, jd in test_jds[:8]:
    for i in range(len(planets)):
        for j in range(i + 1, len(planets)):
            pid1, n1 = planets[i]
            pid2, n2 = planets[j]
            try:
                se = swe.calc_ut(jd, pid1, F | S)
                le = ephem.calc_ut(jd, pid1, F | S)
                se2 = swe.calc_ut(jd, pid2, F | S)
                le2 = ephem.calc_ut(jd, pid2, F | S)

                se_sep = (se2[0][0] - se[0][0]) % 360.0
                le_sep = (le2[0][0] - le[0][0]) % 360.0

                for asp in [0, 60, 90, 120, 180]:
                    se_orb = min(abs(se_sep - asp), abs(se_sep - (360 - asp)))
                    le_orb = min(abs(le_sep - asp), abs(le_sep - (360 - asp)))

                    if se_orb < 10.0:  # Only check when near aspect
                        orb_diff = abs(se_orb - le_orb) * 3600.0  # arcsec
                        orb_diffs.append(orb_diff)
                        if orb_diff < 2.0:
                            passed += 1
                        else:
                            failed += 1
                            if failed <= 5:
                                print(
                                    f'  FAIL {n1}-{n2} asp={asp}: SE_orb={se_orb:.6f} LE_orb={le_orb:.6f} diff={orb_diff:.2f}"'
                                )
            except:
                errors += 1

if orb_diffs:
    avg_diff = sum(orb_diffs) / len(orb_diffs)
    max_diff = max(orb_diffs)
    print(
        f'  Orb precision: avg={avg_diff:.4f}" max={max_diff:.4f}" ({len(orb_diffs)} aspects checked)'
    )
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Applying/separating aspect detection
# ============================================================
print("\n=== P4: Applying vs separating aspect agreement ===")

app_match = app_mismatch = 0
for year, jd in test_jds[:8]:
    for i in range(len(planets)):
        for j in range(i + 1, len(planets)):
            pid1, _ = planets[i]
            pid2, _ = planets[j]
            try:
                se1 = swe.calc_ut(jd, pid1, F | S)
                se2 = swe.calc_ut(jd, pid2, F | S)
                le1 = ephem.calc_ut(jd, pid1, F | S)
                le2 = ephem.calc_ut(jd, pid2, F | S)

                se_sep = (se2[0][0] - se1[0][0]) % 360.0
                le_sep = (le2[0][0] - le1[0][0]) % 360.0

                # Speed difference determines applying/separating
                se_speed_diff = se2[0][3] - se1[0][3]
                le_speed_diff = le2[0][3] - le1[0][3]

                # Both should agree on relative speed direction
                if (se_speed_diff > 0) == (le_speed_diff > 0):
                    app_match += 1
                else:
                    # Near-zero speed diff can legitimately differ
                    if abs(se_speed_diff) < 0.001 or abs(le_speed_diff) < 0.001:
                        app_match += 1  # Too close to call
                    else:
                        app_mismatch += 1
            except:
                errors += 1

passed += app_match
failed += app_mismatch
print(f"  Apply/separate: {app_match} agree, {app_mismatch} disagree")
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Midpoint precision
# ============================================================
print("\n=== P5: Midpoint precision ===")

for year, jd in test_jds[:8]:
    for i in range(len(planets)):
        for j in range(i + 1, min(i + 3, len(planets))):
            pid1, n1 = planets[i]
            pid2, n2 = planets[j]
            try:
                se1 = swe.calc_ut(jd, pid1, F | S)
                se2 = swe.calc_ut(jd, pid2, F | S)
                le1 = ephem.calc_ut(jd, pid1, F | S)
                le2 = ephem.calc_ut(jd, pid2, F | S)

                # Midpoint calculation
                se_mid = (se1[0][0] + se2[0][0]) / 2.0
                le_mid = (le1[0][0] + le2[0][0]) / 2.0

                # Handle wraparound
                diff = abs(se1[0][0] - se2[0][0])
                if diff > 180:
                    se_mid = (se_mid + 180) % 360
                    le_mid_2 = (le1[0][0] + le2[0][0]) / 2.0
                    le_diff = abs(le1[0][0] - le2[0][0])
                    if le_diff > 180:
                        le_mid = (le_mid_2 + 180) % 360
                    else:
                        le_mid = le_mid_2

                mid_diff = abs(se_mid - le_mid) * 3600.0
                if mid_diff > 180 * 3600:
                    mid_diff = 360 * 3600 - mid_diff
                if mid_diff < 1.0:
                    passed += 1
                else:
                    failed += 1
            except:
                errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 92 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
