#!/usr/bin/env python3
"""Round 81: Delta-T at Leap Second Boundaries

Tests Delta-T computation precisely at and around leap second insertion times.
Verifies continuity of Delta-T across leap second boundaries and consistency
with IERS data. Also tests deltat at historical leap seconds (1972-2017).
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

print("=" * 70)
print("ROUND 81: Delta-T at Leap Second Boundaries")
print("=" * 70)

# All leap second dates (UTC) since 1972
# Format: (year, month, day) — leap second inserted at end of this day
leap_seconds = [
    (1972, 6, 30),
    (1972, 12, 31),
    (1973, 12, 31),
    (1974, 12, 31),
    (1975, 12, 31),
    (1976, 12, 31),
    (1977, 12, 31),
    (1978, 12, 31),
    (1979, 12, 31),
    (1981, 6, 30),
    (1982, 6, 30),
    (1983, 6, 30),
    (1985, 6, 30),
    (1987, 12, 31),
    (1989, 12, 31),
    (1990, 12, 31),
    (1992, 6, 30),
    (1993, 6, 30),
    (1994, 6, 30),
    (1995, 12, 31),
    (1997, 6, 30),
    (1998, 12, 31),
    (2005, 12, 31),
    (2008, 12, 31),
    (2012, 6, 30),
    (2015, 6, 30),
    (2016, 12, 31),
]

# ============================================================
# P1: Delta-T at each leap second date
# ============================================================
print("\n=== P1: Delta-T at leap second dates ===")

for year, month, day in leap_seconds:
    jd = swe.julday(year, month, day, 12.0)
    label = f"deltat {year}-{month:02d}-{day:02d}"
    try:
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff_ms = abs(se_dt - le_dt) * 86400000.0  # days to milliseconds
        if diff_ms < 50.0:  # 50ms tolerance
            passed += 1
        else:
            failed += 1
            print(
                f"  FAIL {label}: SE={se_dt:.10f} LE={le_dt:.10f} diff={diff_ms:.1f}ms"
            )
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P2: Delta-T continuity across leap second boundaries
# ============================================================
print("\n=== P2: Delta-T continuity across boundaries ===")

for year, month, day in leap_seconds:
    jd = swe.julday(year, month, day, 23.999)
    # Check day before, day of, day after
    for offset, offset_name in [
        (-1.0, "day_before"),
        (0.0, "day_of"),
        (1.0, "day_after"),
    ]:
        jd_test = jd + offset
        label = f"continuity {year}-{month:02d}-{day:02d} {offset_name}"
        try:
            se_dt = swe.deltat(jd_test)
            le_dt = ephem.deltat(jd_test)
            diff_ms = abs(se_dt - le_dt) * 86400000.0
            if diff_ms < 50.0:
                passed += 1
            else:
                failed += 1
                print(
                    f"  FAIL {label}: SE={se_dt:.10f} LE={le_dt:.10f} diff={diff_ms:.1f}ms"
                )
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P3: Delta-T monotonicity check (should generally increase 1972-2017)
# ============================================================
print("\n=== P3: Delta-T monotonicity (LE internal) ===")

prev_dt = None
monotonic_pass = 0
monotonic_fail = 0
for year in range(1972, 2018):
    jd = swe.julday(year, 7, 1, 12.0)
    le_dt = ephem.deltat(jd)
    if prev_dt is not None:
        if le_dt >= prev_dt - 0.000001:  # Allow tiny rounding
            monotonic_pass += 1
        else:
            monotonic_fail += 1
            print(
                f"  FAIL monotonicity at {year}: dt={le_dt * 86400:.3f}s < prev={prev_dt * 86400:.3f}s"
            )
    prev_dt = le_dt

passed += monotonic_pass
failed += monotonic_fail
print(f"  Monotonicity: {monotonic_pass} passed, {monotonic_fail} failed")
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P4: Delta-T at specific known values (IERS Bulletin A)
# ============================================================
print("\n=== P4: Delta-T at known IERS values ===")

# Known Delta-T values (seconds) from IERS
known_dt = [
    (1972, 1, 1, 42.23),
    (1975, 1, 1, 45.48),
    (1980, 1, 1, 50.54),
    (1985, 1, 1, 54.34),
    (1990, 1, 1, 56.86),
    (1995, 1, 1, 60.78),
    (2000, 1, 1, 63.83),
    (2005, 1, 1, 64.69),
    (2010, 1, 1, 66.07),
    (2015, 1, 1, 67.64),
    (2017, 1, 1, 68.59),
]

for year, month, day, expected_dt_s in known_dt:
    jd = swe.julday(year, month, day, 0.0)
    label = f"IERS {year}-{month:02d}-{day:02d} (expected {expected_dt_s:.2f}s)"
    try:
        le_dt_s = ephem.deltat(jd) * 86400.0
        diff_s = abs(le_dt_s - expected_dt_s)
        if diff_s < 1.0:  # Within 1 second of IERS value
            passed += 1
        else:
            failed += 1
            print(f"  FAIL {label}: LE={le_dt_s:.3f}s diff={diff_s:.3f}s")
    except Exception as e:
        errors += 1
        print(f"  ERROR {label}: {e}")

print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P5: Delta-T SE vs LE fine grid (every month 1970-2025)
# ============================================================
print("\n=== P5: Delta-T fine grid monthly 1970-2025 ===")

fine_pass = 0
fine_fail = 0
for year in range(1970, 2026):
    for month in range(1, 13):
        jd = swe.julday(year, month, 15, 12.0)
        try:
            se_dt = swe.deltat(jd)
            le_dt = ephem.deltat(jd)
            diff_ms = abs(se_dt - le_dt) * 86400000.0
            # Modern era should be very close
            if diff_ms < 100.0:  # 100ms tolerance
                fine_pass += 1
            else:
                fine_fail += 1
                if fine_fail <= 10:
                    print(
                        f"  FAIL {year}-{month:02d}: SE={se_dt * 86400:.3f}s LE={le_dt * 86400:.3f}s diff={diff_ms:.1f}ms"
                    )
        except Exception as e:
            errors += 1

passed += fine_pass
failed += fine_fail
if fine_fail > 10:
    print(f"  ... and {fine_fail - 10} more failures")
print(f"  Fine grid: {fine_pass} passed, {fine_fail} failed")
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P6: Delta-T at historical dates (pre-1972)
# ============================================================
print("\n=== P6: Delta-T historical (1800-1972) ===")

hist_pass = 0
hist_fail = 0
for year in range(1800, 1972, 10):
    jd = swe.julday(year, 1, 1, 12.0)
    try:
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff_s = abs(se_dt - le_dt) * 86400.0
        # Historical Delta-T can differ more between models
        if diff_s < 20.0:  # 20 second tolerance for pre-1972
            hist_pass += 1
        else:
            hist_fail += 1
            print(
                f"  FAIL {year}: SE={se_dt * 86400:.2f}s LE={le_dt * 86400:.2f}s diff={diff_s:.2f}s"
            )
    except Exception as e:
        errors += 1

passed += hist_pass
failed += hist_fail
print(f"  Historical: {hist_pass} passed, {hist_fail} failed")
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P7: Delta-T future prediction (2025-2050)
# ============================================================
print("\n=== P7: Delta-T future prediction (2025-2050) ===")

fut_pass = 0
fut_fail = 0
for year in range(2025, 2051):
    jd = swe.julday(year, 7, 1, 12.0)
    try:
        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)
        diff_s = abs(se_dt - le_dt) * 86400.0
        # Future predictions can diverge — allow up to 30s
        if diff_s < 30.0:
            fut_pass += 1
        else:
            fut_fail += 1
            print(
                f"  FAIL {year}: SE={se_dt * 86400:.2f}s LE={le_dt * 86400:.2f}s diff={diff_s:.2f}s"
            )
    except Exception as e:
        errors += 1

passed += fut_pass
failed += fut_fail
print(f"  Future: {fut_pass} passed, {fut_fail} failed")
print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# P8: Impact on planetary positions at leap second boundaries
# ============================================================
print("\n=== P8: Position impact at leap second boundaries ===")

for year, month, day in leap_seconds[-5:]:  # Last 5 leap seconds
    jd = swe.julday(year, month, day, 12.0)
    for body, name in [(0, "Sun"), (1, "Moon"), (4, "Mars")]:
        label = f"{name} {year}-{month:02d}-{day:02d}"
        try:
            se_pos = swe.calc_ut(jd, body, swe.FLG_SWIEPH | swe.FLG_SPEED)
            le_pos = ephem.calc_ut(jd, body, 2 | 256)
            diff_arcsec = abs(se_pos[0][0] - le_pos[0]) * 3600.0
            if diff_arcsec > 180 * 3600:
                diff_arcsec = 360 * 3600 - diff_arcsec
            if diff_arcsec < 1.0:
                passed += 1
            else:
                failed += 1
                print(
                    f'  FAIL {label}: SE={se_pos[0][0]:.6f} LE={le_pos[0]:.6f} diff={diff_arcsec:.2f}"'
                )
        except Exception as e:
            errors += 1
            print(f"  ERROR {label}: {e}")

print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

# ============================================================
# Summary
# ============================================================
total = passed + failed
print("\n" + "=" * 70)
print(f"ROUND 81 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print("=" * 70)
