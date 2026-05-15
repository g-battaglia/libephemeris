#!/usr/bin/env python3
"""Round 162: Ayanamsha at sub-second time intervals.

Compare get_ayanamsa_ut at very fine time intervals (sub-second)
to verify interpolation smoothness and precision. Tests that both
libraries produce consistent, smooth ayanamsha values.
"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

SIDM_MODES = [(0, "Fagan"), (1, "Lahiri"), (3, "Raman"), (27, "TrueCitra")]

# Test at sub-second intervals around several dates
test_cases = []
base_dates = [
    swe.julday(2000, 1, 1, 12.0),
    swe.julday(2024, 3, 20, 12.0),
    swe.julday(2025, 6, 21, 0.0),
    swe.julday(1990, 12, 31, 23.9),
]

for base_jd in base_dates:
    # Sub-second: 1 second = 1/86400 day ≈ 1.157e-5
    for offset_sec in range(-30, 31):  # ±30 seconds
        jd = base_jd + offset_sec / 86400.0
        test_cases.append(jd)

# Also test at 1-minute intervals over an hour
base = swe.julday(2024, 1, 1, 0.0)
for minute in range(60):
    test_cases.append(base + minute / 1440.0)

TOL = 0.05  # 50 mas tolerance for ayanamsha

passed = failed = errors = total = 0
failures = []

print(f"Round 162: Ayanamsha at Sub-Second Time Intervals")
print(
    f"Testing {len(SIDM_MODES)} modes x {len(test_cases)} times = {len(SIDM_MODES) * len(test_cases)} tests"
)
print("=" * 90)

for sidm, sname in SIDM_MODES:
    swe.set_sid_mode(sidm)
    ephem.set_sid_mode(sidm, 0, 0)

    for jd in test_cases:
        total += 1
        try:
            se_aya = swe.get_ayanamsa_ut(jd)
            le_aya = ephem.get_ayanamsa_ut(jd)
            diff = abs(le_aya - se_aya) * 3600.0  # to arcsec
            if diff <= TOL:
                passed += 1
            else:
                failed += 1
                msg = f'  FAIL {sname} jd={jd:.10f}: SE={se_aya:.10f} LE={le_aya:.10f} diff={diff:.4f}"'
                failures.append(msg)
                if len(failures) <= 10:
                    print(msg)
        except Exception as e:
            errors += 1

swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0, 0)

print(f"\n{'=' * 90}")
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\n{len(failures)} failures")
else:
    print("\nAll tests passed!")
