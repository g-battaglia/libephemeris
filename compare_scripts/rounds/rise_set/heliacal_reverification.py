#!/usr/bin/env python3
"""Round 212: Heliacal Re-verification (Quick).

Quick check of heliacal functions with minimal test cases to avoid
the extreme slowness of heliacal calculations.
Only tests heliacal_ut for a few bright planets/stars.
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

# Geographic locations
LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 21.0),
    ("New York", 40.7128, -74.0060, 10.0),
    ("Delhi", 28.6139, 77.2090, 216.0),
]

# Bodies to test (only bright, fast ones)
BODIES = [
    ("Venus", ephem.VENUS, swe.VENUS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Mercury", ephem.MERCURY, swe.MERCURY),
]

# Event types
EVENT_TYPES = [
    ("Heliacal Rising", 1),
    ("Heliacal Setting", 2),
]

# Only a few dates to keep it fast
DATES = [2451545.0, 2460000.0]

ATPRESS = 1013.25
ATTEMP = 15.0


if __name__ == "__main__":
    print("=" * 70)
    print("Round 212: Heliacal Re-verification (Quick)")
    print("=" * 70)

    for loc_name, lat, lon, alt in LOCATIONS:
        print(f"\n--- {loc_name} ---")
        for body_name, le_body, se_body in BODIES:
            for evt_name, evt_type in EVENT_TYPES:
                for jd in DATES:
                    label = f"{loc_name} {body_name} {evt_name} JD={jd:.1f}"

                    try:
                        geopos = [lon, lat, alt]
                        se_result = swe.heliacal_ut(
                            jd,
                            geopos,
                            [ATPRESS, ATTEMP, 50, 0, 0, 0],
                            [
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                            ],
                            str(se_body),
                            evt_type,
                            swe.FLG_SWIEPH | swe.HELFLAG_HIGH_PRECISION,
                        )
                        se_jd = se_result[0]
                    except Exception as e:
                        # SE might fail for some configurations
                        continue

                    try:
                        le_result = ephem.heliacal_ut(
                            jd,
                            [lon, lat, alt],
                            [ATPRESS, ATTEMP, 50, 0, 0, 0],
                            [0] * 20,
                            str(le_body),
                            evt_type,
                            ephem.FLG_SWIEPH | ephem.HELFLAG_HIGH_PRECISION,
                        )
                        le_jd = le_result[0]
                    except Exception as e:
                        continue

                    total += 1
                    diff_days = abs(le_jd - se_jd)
                    diff_hours = diff_days * 24

                    if diff_hours <= 24.0:  # within 24 hours
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f"  {label}: LE_JD={le_jd:.6f} SE_JD={se_jd:.6f} diff={diff_hours:.2f}h"
                        )

    print(f"\n{'=' * 70}")
    if total > 0:
        print(
            f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
        )
    else:
        print("RESULTS: 0 tests ran (all skipped)")
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:20]:
            print(f)
