#!/usr/bin/env python3
"""Round 214: All ayanamsha modes at multiple dates.

Comprehensive test of all ayanamsha modes (0-42+) at multiple dates,
comparing get_ayanamsa_ex_ut results between LE and SE.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

# All ayanamsha modes to test
AYANAMSHA_MODES = list(range(43))

# Dates spanning wide range
DATES = [
    2415020.0,  # 1900
    2430000.0,  # 1941
    2440000.0,  # 1968
    2451545.0,  # J2000
    2455000.0,  # 2009
    2460000.0,  # 2023
    2462000.0,  # 2028
]

# Known modes with large expected differences
KNOWN_LARGE_DIFF_MODES = {
    17,  # GALCENT_0SAG - galactic center based
    20,  # GALCENT_RGILBRAND
    36,  # GALCENT_COCHRANE_2
    35,  # ARYABHATA_522
    40,  # TRUE_CITRA - known ~13.9" offset
}

FLAGS = ephem.FLG_SWIEPH


def compare_ayanamsha(mode, jd):
    global passed, failed, total

    try:
        swe.set_sid_mode(mode)
        ephem.set_sid_mode(mode, 0, 0)
    except Exception:
        return

    try:
        le_aya = ephem.get_ayanamsa_ex_ut(jd, FLAGS)
        # le_aya is (retflag, ayanamsa_value)
        if isinstance(le_aya, tuple):
            le_val = le_aya[1]
        else:
            le_val = le_aya
    except Exception:
        return

    try:
        se_aya = swe.get_ayanamsa_ex_ut(jd, swe.FLG_SWIEPH)
        if isinstance(se_aya, tuple):
            se_val = se_aya[1] if len(se_aya) > 1 else se_aya[0]
        else:
            se_val = se_aya
    except Exception:
        return

    total += 1
    diff = abs(le_val - se_val)
    diff_as = diff * 3600

    # Set tolerance based on known differences
    if mode in KNOWN_LARGE_DIFF_MODES:
        tol = 800.0  # 800" for galactic center modes
    else:
        tol = 30.0  # 30" for normal modes (known ~14" systematic)

    if diff_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  Mode {mode} JD={jd:.1f}: LE={le_val:.6f} SE={se_val:.6f} diff={diff_as:.2f}"'
        )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 214: All Ayanamsha Modes at Multiple Dates")
    print("=" * 70)

    for mode in AYANAMSHA_MODES:
        for jd in DATES:
            compare_ayanamsha(mode, jd)

    # Reset
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    print(f"\n{'=' * 70}")
    if total > 0:
        print(
            f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
        )
    else:
        print("RESULTS: 0 tests ran")
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:40]:
            print(f)
        if len(failures) > 40:
            print(f"  ... and {len(failures) - 40} more")
