#!/usr/bin/env python3
"""Round 155: Aberration correction magnitude verification.

Compare positions with and without FLG_NOABERR to verify that the
aberration correction magnitude (~20.5" for annual aberration) is
consistent between libephemeris and pyswisseph.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_NOABERR = 1024

BODIES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

test_dates = []
for year in [1950, 1980, 2000, 2010, 2020, 2024, 2025, 2050]:
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 1, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

TOL_ABERR_DIFF = 5.0  # arcsec tolerance on aberration difference itself

passed = failed = errors = total = 0
failures = []

print(f"Round 155: Aberration Correction Magnitude Verification")
print(f"Testing {len(BODIES)} bodies x {len(test_dates)} dates")
print("=" * 90)

for label, jd in test_dates:
    for body, bname in BODIES.items():
        try:
            # With aberration (default)
            se_with = swe.calc_ut(jd, body, FLG_SPEED)[0]
            le_with = ephem.calc_ut(jd, body, FLG_SPEED)[0]

            # Without aberration
            se_without = swe.calc_ut(jd, body, FLG_SPEED | FLG_NOABERR)[0]
            le_without = ephem.calc_ut(jd, body, FLG_SPEED | FLG_NOABERR)[0]

            # Aberration correction = position_with - position_without
            se_aberr_lon = (se_with[0] - se_without[0]) * 3600.0
            le_aberr_lon = (le_with[0] - le_without[0]) * 3600.0
            se_aberr_lat = (se_with[1] - se_without[1]) * 3600.0
            le_aberr_lat = (le_with[1] - le_without[1]) * 3600.0

            # Compare aberration magnitudes
            for comp, se_a, le_a in [
                ("aberr_lon", se_aberr_lon, le_aberr_lon),
                ("aberr_lat", se_aberr_lat, le_aberr_lat),
            ]:
                total += 1
                diff = abs(le_a - se_a)
                if diff <= TOL_ABERR_DIFF:
                    passed += 1
                else:
                    failed += 1
                    msg = f'  FAIL {label} {bname} {comp}: SE={se_a:.4f}" LE={le_a:.4f}" diff={diff:.4f}"'
                    failures.append(msg)
                    if len(failures) <= 20:
                        print(msg)

            # Also compare the NOABERR positions directly
            for i, cname in enumerate(["lon", "lat"]):
                total += 1
                diff = abs(le_without[i] - se_without[i]) * 3600.0
                tol = 1.0  # 1" for NOABERR positions
                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    msg = f'  FAIL {label} {bname} NOABERR_{cname}: SE={se_without[i]:.8f} LE={le_without[i]:.8f} diff={diff:.4f}"'
                    failures.append(msg)
                    if len(failures) <= 20:
                        print(msg)

        except Exception as e:
            errors += 1

print()
print("=" * 90)
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\n{len(failures)} failures")
else:
    print("\nAll tests passed!")
