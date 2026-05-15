#!/usr/bin/env python3
"""Round 83: NONUT/NOABERR/TRUEPOS Flag Deep Sweep

Deep verification of special calculation flags that modify the position pipeline:
- FLG_NONUT (64): No nutation
- FLG_NOABERR (1024): No aberration correction
- FLG_TRUEPOS (16): True geometric position (no light-time, no aberration)
- Combined flags
- All planets across 50-year sweep
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0

FLG_SWIEPH = 2
FLG_SPEED = 256
FLG_TRUEPOS = 16
FLG_J2000 = 32
FLG_NONUT = 64
FLG_NOABERR = 1024
FLG_EQUATORIAL = 2048


def check_pos(label, se_pos, le_pos, lon_tol=1.0, lat_tol=1.0, dist_tol=0.001):
    """Compare position tuples. Tolerances in arcsec for lon/lat, AU for dist."""
    global passed, failed
    for i, (name, tol) in enumerate(
        [("lon", lon_tol), ("lat", lat_tol), ("dist", dist_tol)]
    ):
        if i < 2:
            diff = abs(se_pos[i] - le_pos[i]) * 3600.0
            if diff > 180 * 3600:
                diff = 360 * 3600 - diff
            if diff < tol:
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label} {name}: SE={se_pos[i]:.6f} LE={le_pos[i]:.6f} diff={diff:.2f}"'
                )
        else:
            if se_pos[i] == 0 and le_pos[i] == 0:
                passed += 1
                continue
            ratio = se_pos[i] / le_pos[i] if le_pos[i] != 0 else 999
            if abs(ratio - 1.0) < tol:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL {label} {name}: SE={se_pos[i]:.8f} LE={le_pos[i]:.8f} ratio={ratio:.8f}"
                )


print("=" * 70)
print("ROUND 83: NONUT/NOABERR/TRUEPOS Flag Deep Sweep")
print("=" * 70)

bodies = [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MERCURY, "Mercury"),
    (swe.VENUS, "Venus"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
    (swe.URANUS, "Uranus"),
    (swe.NEPTUNE, "Neptune"),
    (swe.PLUTO, "Pluto"),
]

# Test dates: 25 years, every 2 years
test_jds = []
for year in range(1980, 2030, 2):
    jd = swe.julday(year, 1, 15, 12.0)
    test_jds.append((year, jd))

# ============================================================
# P1: NONUT flag (ecliptic of date without nutation = mean ecliptic)
# ============================================================
print("\n=== P1: NONUT flag all planets ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_NONUT

for body_id, name in bodies:
    body_pass = 0
    body_fail = 0
    for year, jd in test_jds:
        label = f"{name} {year} NONUT"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            check_pos(label, se_pos, le_pos)
        except Exception as e:
            errors += 1
            if errors <= 5:
                print(f"  ERROR {label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: NOABERR flag
# ============================================================
print("\n=== P2: NOABERR flag all planets ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_NOABERR

for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year} NOABERR"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            check_pos(label, se_pos, le_pos)
        except Exception as e:
            errors += 1

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: TRUEPOS flag
# ============================================================
print("\n=== P3: TRUEPOS flag all planets ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS

for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year} TRUEPOS"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            # TRUEPOS can have larger diff due to light-time treatment
            check_pos(label, se_pos, le_pos, lon_tol=2.0, lat_tol=2.0)
        except Exception as e:
            errors += 1

print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: NONUT + NOABERR combined
# ============================================================
print("\n=== P4: NONUT + NOABERR combined ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_NONUT | FLG_NOABERR

for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year} NONUT+NOABERR"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            check_pos(label, se_pos, le_pos)
        except Exception as e:
            errors += 1

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: TRUEPOS + NONUT combined
# ============================================================
print("\n=== P5: TRUEPOS + NONUT combined ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS | FLG_NONUT

for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year} TRUEPOS+NONUT"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            check_pos(label, se_pos, le_pos, lon_tol=2.0, lat_tol=2.0)
        except Exception as e:
            errors += 1

print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: J2000 + NONUT (should be identical since J2000 has no nutation)
# ============================================================
print("\n=== P6: J2000 + NONUT ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_J2000 | FLG_NONUT

for body_id, name in bodies:
    for year, jd in test_jds:
        label = f"{name} {year} J2000+NONUT"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            check_pos(label, se_pos, le_pos)
        except Exception as e:
            errors += 1

print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: EQUATORIAL + NONUT (mean equator)
# ============================================================
print("\n=== P7: EQUATORIAL + NONUT (mean equator) ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL | FLG_NONUT

for body_id, name in [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
    (swe.SATURN, "Saturn"),
]:
    for year, jd in test_jds:
        label = f"{name} {year} EQ+NONUT"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            # RA in degrees, Dec in degrees
            check_pos(label, se_pos, le_pos, lon_tol=1.5, lat_tol=1.5)
        except Exception as e:
            errors += 1

print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P8: EQUATORIAL + NOABERR
# ============================================================
print("\n=== P8: EQUATORIAL + NOABERR ===")

base_flags = FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL | FLG_NOABERR

for body_id, name in [
    (swe.SUN, "Sun"),
    (swe.MOON, "Moon"),
    (swe.MARS, "Mars"),
    (swe.JUPITER, "Jupiter"),
]:
    for year, jd in test_jds:
        label = f"{name} {year} EQ+NOABERR"
        try:
            se = swe.calc_ut(jd, body_id, base_flags)
            le = ephem.calc_ut(jd, body_id, base_flags)
            se_pos = se[0][:3]
            le_pos = le[0][:3]
            check_pos(label, se_pos, le_pos, lon_tol=1.5, lat_tol=1.5)
        except Exception as e:
            errors += 1

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P9: Verify NONUT effect magnitude (~9" nutation amplitude)
# ============================================================
print("\n=== P9: Nutation effect magnitude verification ===")

for year, jd in test_jds[:5]:
    for body_id, name in [(swe.SUN, "Sun"), (swe.MARS, "Mars")]:
        label = f"{name} {year} nutation_effect"
        try:
            # With nutation
            le_nut = ephem.calc_ut(jd, body_id, FLG_SWIEPH | FLG_SPEED)
            # Without nutation
            le_nonut = ephem.calc_ut(
                jd, body_id, FLG_SWIEPH | FLG_SPEED | FLG_NONUT
            )

            diff_arcsec = abs(le_nut[0][0] - le_nonut[0][0]) * 3600.0
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec

            # Nutation in longitude typically 0-20" (max ~17.2")
            if diff_arcsec < 25.0:
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: nutation effect={diff_arcsec:.2f}" (expected < 25")'
                )
        except Exception as e:
            errors += 1

print(f"  After P9: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P10: Verify NOABERR effect magnitude (~20" for Sun)
# ============================================================
print("\n=== P10: Aberration effect magnitude verification ===")

for year, jd in test_jds[:5]:
    for body_id, name in [(swe.SUN, "Sun"), (swe.MARS, "Mars")]:
        label = f"{name} {year} aberration_effect"
        try:
            # With aberration
            le_aberr = ephem.calc_ut(jd, body_id, FLG_SWIEPH | FLG_SPEED)
            # Without aberration
            le_noaberr = ephem.calc_ut(
                jd, body_id, FLG_SWIEPH | FLG_SPEED | FLG_NOABERR
            )

            diff_arcsec = abs(le_aberr[0][0] - le_noaberr[0][0]) * 3600.0
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec

            # Annual aberration ~20.5" max for Sun, varies for planets
            if diff_arcsec < 40.0:  # generous for planets
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: aberration effect={diff_arcsec:.2f}" (expected < 40")'
                )
        except Exception as e:
            errors += 1

print(f"  After P10: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P11: Speed values with flags
# ============================================================
print("\n=== P11: Speed consistency with special flags ===")

for body_id, name in [(swe.SUN, "Sun"), (swe.MOON, "Moon"), (swe.MARS, "Mars")]:
    for year, jd in test_jds[:10]:
        for flag_extra, flag_name in [
            (FLG_NONUT, "NONUT"),
            (FLG_NOABERR, "NOABERR"),
            (FLG_TRUEPOS, "TRUEPOS"),
        ]:
            flags = FLG_SWIEPH | FLG_SPEED | flag_extra
            label = f"{name} {year} speed {flag_name}"
            try:
                se = swe.calc_ut(jd, body_id, flags)
                le = ephem.calc_ut(jd, body_id, flags)
                se_speed = se[0][3]  # lon speed
                le_speed = le[0][3]
                diff_speed = abs(se_speed - le_speed) * 3600.0  # arcsec/day
                if diff_speed < 1.0:  # 1"/day tolerance
                    passed += 1
                else:
                    failed += 1
                    print(
                        f'  FAIL {label}: SE={se_speed:.6f} LE={le_speed:.6f} diff={diff_speed:.2f}"/day'
                    )
            except Exception as e:
                errors += 1

print(f"  After P11: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# Summary
# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 83 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
