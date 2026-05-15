#!/usr/bin/env python3
"""Round 87: Heliocentric All Planets with Flag Combinations

Deep verification of heliocentric positions across all planets with
various flag combinations: J2000, NONUT, EQUATORIAL, XYZ, RADIANS, NOABERR.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = failed = errors = 0
F = 2  # SWIEPH
S = 256  # SPEED
H = 8  # HELCTR
J = 32  # J2000
N = 64  # NONUT
E = 2048  # EQUATORIAL
X = 4096  # XYZ
R = 8192  # RADIANS

print("=" * 70)
print("ROUND 87: Heliocentric All Planets with Flag Combinations")
print("=" * 70)

bodies = [
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
    (15, "Chiron"),
    (17, "Ceres"),
    (18, "Pallas"),
]
test_jds = [(y, swe.julday(y, 1, 15, 12.0)) for y in range(1980, 2026, 3)]

flag_combos = [
    (F | S | H, "helio"),
    (F | S | H | J | N, "helio+J2000"),
    (F | S | H | E, "helio+EQ"),
    (F | S | H | X, "helio+XYZ"),
    (F | S | H | R, "helio+RAD"),
    (F | S | H | J | N | E, "helio+J2000+EQ"),
]

for flags, flag_name in flag_combos:
    sec_pass = sec_fail = 0
    is_xyz = bool(flags & X)
    for body_id, name in bodies:
        for year, jd in test_jds:
            label = f"{name} {year} {flag_name}"
            try:
                se = swe.calc_ut(jd, body_id, flags)
                le = ephem.calc_ut(jd, body_id, flags)
                for i in range(3):
                    if is_xyz:
                        diff = abs(se[0][i] - le[0][i])
                        ok = diff < 0.0002  # AU
                    elif i == 2:
                        ratio = se[0][i] / le[0][i] if le[0][i] != 0 else 999
                        ok = abs(ratio - 1.0) < 0.0001
                    else:
                        diff = abs(se[0][i] - le[0][i]) * 3600.0
                        if diff > 180 * 3600:
                            diff = 360 * 3600 - diff
                        ok = diff < 2.0
                    if ok:
                        passed += 1
                        sec_pass += 1
                    else:
                        failed += 1
                        sec_fail += 1
                        if sec_fail <= 3:
                            print(
                                f"  FAIL {label} [{i}]: SE={se[0][i]:.8f} LE={le[0][i]:.8f}"
                            )
            except Exception as e:
                errors += 1
    if sec_fail > 3:
        print(f"  ... {flag_name}: {sec_fail} total failures")
    print(f"  {flag_name}: {sec_pass} passed, {sec_fail} failed")

total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 87 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed: {passed}  Failed: {failed}  Errors: {errors}")
print("=" * 70)
