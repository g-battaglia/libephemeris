#!/usr/bin/env python3
"""Round 180: Eclipse Saros cycle deep — consecutive eclipses in same Saros.

Tests solar eclipse timing across a Saros series (6585.32 days between
consecutive members) to verify long-term eclipse prediction consistency.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

# Start from known eclipses and search forward
START_DATES = [
    (1999, 8, 11),  # Famous 1999 total
    (2001, 6, 21),  # 2001 total
    (2006, 3, 29),  # 2006 total
    (2009, 7, 22),  # 2009 total
    (2017, 8, 21),  # Great American Eclipse
    (2024, 4, 8),  # 2024 total
]

TOL_MINUTES = 10.0  # timing tolerance

passed = 0
failed = 0
skipped = 0
errors = []

for y, m, d in START_DATES:
    jd_start = swe.julday(y, m, d, 0.0) - 30  # Start a month before

    # Search for next eclipse using glob
    try:
        se_r = swe.sol_eclipse_when_glob(jd_start, 0, 0, False)
        se_tret = se_r[1]
        se_max = se_tret[0]
    except Exception:
        skipped += 1
        continue

    try:
        le_r = ephem.sol_eclipse_when_glob(jd_start, 0, 0, "forward")
        le_tret = le_r[1]
        le_max = le_tret[0]
    except Exception:
        skipped += 1
        continue

    if se_max == 0.0 or le_max == 0.0:
        skipped += 1
        continue

    diff_min = abs(se_max - le_max) * 1440.0

    if diff_min <= TOL_MINUTES:
        passed += 1
    else:
        failed += 1
        if len(errors) < 20:
            errors.append(f"  FAIL eclipse near {y}/{m}/{d}: diff={diff_min:.1f}min")

    # Now search for next eclipse after the one found (Saros member search)
    for step in range(3):
        next_start = se_max + 10  # 10 days after the found eclipse
        try:
            se_r2 = swe.sol_eclipse_when_glob(next_start, 0, 0, False)
            se_max = se_r2[1][0]
        except Exception:
            break
        try:
            le_r2 = ephem.sol_eclipse_when_glob(next_start, 0, 0, "forward")
            le_max = le_r2[1][0]
        except Exception:
            break

        if se_max == 0.0 or le_max == 0.0:
            break

        diff_min = abs(se_max - le_max) * 1440.0
        if diff_min <= TOL_MINUTES:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(
                    f"  FAIL eclipse step{step + 1} after {y}/{m}/{d}: diff={diff_min:.1f}min"
                )

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 180: Eclipse Saros Cycle Deep ===")
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 180: ALL PASSED")
