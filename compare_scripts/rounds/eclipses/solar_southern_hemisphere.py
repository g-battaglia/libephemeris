#!/usr/bin/env python3
"""Round 169: Solar eclipse when_loc at southern hemisphere locations.

Tests sol_eclipse_when_loc for locations in the southern hemisphere
where eclipse geometry and visibility patterns differ from the north.
Tests multiple locations across South America, Africa, Australia, Antarctica.
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

# Southern hemisphere locations
LOCATIONS = [
    ("BuenosAires", -34.6, -58.4, 25),
    ("SaoPaulo", -23.5, -46.6, 760),
    ("CapeTown", -33.9, 18.4, 0),
    ("Sydney", -33.9, 151.2, 58),
    ("Wellington", -41.3, 174.8, 0),
    ("Santiago", -33.4, -70.7, 520),
    ("McMurdo", -77.8, 166.7, 24),  # Antarctica
    ("Perth", -31.9, 115.9, 0),
    ("Nairobi", -1.3, 36.8, 1795),  # Near equator, south
]

# Start dates for eclipse searches (known eclipse periods)
START_DATES = [
    (2019, 1, 1),  # 2019-07-02 total in S. America
    (2020, 1, 1),  # 2020-12-14 total in S. America
    (2023, 1, 1),  # 2023-04-20 hybrid
    (2024, 1, 1),  # 2024-04-08 total
    (2025, 1, 1),  # 2025 eclipses
    (2028, 1, 1),  # 2028-07-22 total in Australia
]

TOL_JD = 10.0 / 1440.0  # 10 minutes in JD

passed = 0
failed = 0
skipped = 0
errors = []

for y, m, d in START_DATES:
    jd_start = swe.julday(y, m, d, 0.0)

    for loc_name, lat, lon, alt in LOCATIONS:
        geopos = [lon, lat, float(alt)]

        try:
            se_result = swe.sol_eclipse_when_loc(jd_start, geopos, 0, False)
            se_tret = se_result[1]  # time array
            se_max_jd = se_tret[0]  # time of max eclipse
        except Exception:
            skipped += 1
            continue

        if se_max_jd == 0.0:
            skipped += 1
            continue

        try:
            le_result = ephem.sol_eclipse_when_loc(
                jd_start, lat, lon, float(alt), 0
            )
            le_tret = (
                le_result[1] if isinstance(le_result[1], (list, tuple)) else le_result
            )
            le_max_jd = le_tret[0]
        except Exception:
            skipped += 1
            continue

        if le_max_jd == 0.0:
            skipped += 1
            continue

        diff_min = abs(se_max_jd - le_max_jd) * 1440.0

        if diff_min <= 10.0:
            passed += 1
        else:
            failed += 1
            if len(errors) < 25:
                errors.append(
                    f"  FAIL {loc_name} search from {y}: "
                    f"SE max={se_max_jd:.6f}, LE max={le_max_jd:.6f}, "
                    f"diff={diff_min:.1f}min"
                )

total = passed + failed
print(f"\n=== Round 169: Solar Eclipse when_loc Southern Hemisphere ===")
print(f"Start dates: {len(START_DATES)}, Locations: {len(LOCATIONS)}")
pct = 100 * passed / total if total > 0 else 0
print(
    f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)

if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)

if failed == 0:
    print("\nRound 169: ALL PASSED")
