#!/usr/bin/env python3
"""Round 151: Planetary cazimi/combust boundary detection.

Test planet positions near conjunction with Sun to verify sub-degree precision
at cazimi (0°17') and combust (varies by planet) boundaries. Tests that
both libraries agree on exact elongation from Sun.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
SUN = 0

BODIES = {2: "Mercury", 3: "Venus", 4: "Mars", 5: "Jupiter", 6: "Saturn"}

# Approximate inferior/superior conjunction dates for inner planets
# and conjunction dates for outer planets (2020-2025)
# We'll sample densely around these dates to catch near-cazimi moments
test_cases = []

# Dense sampling: every day for 2 months around known conjunction dates
conj_seeds = [
    # Mercury inferior conjunctions (approx)
    (2020, 2, 10),
    (2020, 6, 1),
    (2020, 10, 25),
    (2021, 1, 23),
    (2021, 5, 14),
    (2021, 10, 9),
    (2022, 1, 7),
    (2022, 5, 1),
    (2022, 9, 23),
    (2023, 3, 17),
    (2023, 8, 1),
    (2023, 12, 22),
    (2024, 3, 1),
    (2024, 7, 15),
    (2024, 12, 6),
    # Venus conjunction
    (2022, 1, 9),
    (2023, 8, 13),
    (2024, 6, 4),
    # Mars conjunction
    (2021, 10, 8),
    (2023, 11, 18),
    # Jupiter conjunction
    (2021, 1, 29),
    (2022, 3, 5),
    (2023, 4, 11),
    (2024, 12, 7),
    # Saturn conjunction
    (2021, 1, 24),
    (2022, 2, 4),
    (2023, 2, 16),
    (2024, 2, 28),
]

for year, month, day in conj_seeds:
    # Sample ±15 days around conjunction
    base_jd = swe.julday(year, month, day, 12.0)
    for offset in range(-15, 16):
        jd = base_jd + offset
        test_cases.append((f"{year}-{month:02d}-{day:02d}+{offset:+d}d", jd))

# Tolerances
TOL_LON = 1.0  # arcsec
TOL_LAT = 1.0  # arcsec
TOL_DIST = 1e-6  # AU
TOL_SPD = 1.0  # "/day

passed = failed = errors = total = 0
failures = []

print(f"Round 151: Planetary Cazimi/Combust Boundary Detection")
print(
    f"Testing {len(BODIES)} planets x {len(test_cases)} dates = {len(BODIES) * len(test_cases)} combos"
)
print("=" * 90)

for label, jd in test_cases:
    for body, bname in BODIES.items():
        try:
            # Get Sun position from both
            se_sun = swe.calc_ut(jd, SUN, FLG_SPEED)[0]
            le_sun = ephem.calc_ut(jd, SUN, FLG_SPEED)[0]

            # Get planet position from both
            se_pl = swe.calc_ut(jd, body, FLG_SPEED)[0]
            le_pl = ephem.calc_ut(jd, body, FLG_SPEED)[0]

            # Compare planet positions
            for i, (cname, tol, unit) in enumerate(
                [
                    ("lon", TOL_LON, '"'),
                    ("lat", TOL_LAT, '"'),
                    ("dist", TOL_DIST, "AU"),
                    ("lon_spd", TOL_SPD, '"/d'),
                    ("lat_spd", TOL_SPD, '"/d'),
                    ("dist_spd", 1e-6, "AU/d"),
                ]
            ):
                total += 1
                se_val = se_pl[i]
                le_val = le_pl[i]

                if i <= 1 or i == 3 or i == 4:
                    diff = abs(le_val - se_val) * 3600.0
                else:
                    diff = abs(le_val - se_val)

                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    msg = f"  FAIL {label} {bname} {cname}: SE={se_val:.8f} LE={le_val:.8f} diff={diff:.6f}{unit}"
                    failures.append(msg)
                    if len(failures) <= 15:
                        print(msg)

            # Also compare elongation (separation from Sun)
            total += 1
            se_elong = abs(se_pl[0] - se_sun[0])
            if se_elong > 180:
                se_elong = 360 - se_elong
            le_elong = abs(le_pl[0] - le_sun[0])
            if le_elong > 180:
                le_elong = 360 - le_elong

            diff_elong = abs(le_elong - se_elong) * 3600.0
            if diff_elong <= 2.0:  # 2" tolerance for derived elongation
                passed += 1
            else:
                failed += 1
                msg = f'  FAIL {label} {bname} elongation: SE={se_elong:.6f}° LE={le_elong:.6f}° diff={diff_elong:.2f}"'
                failures.append(msg)
                if len(failures) <= 15:
                    print(msg)

        except Exception as e:
            errors += 1

print()
print("=" * 90)
print(
    f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
)
if failures:
    print(f"\nTotal failures: {len(failures)}")
    cats = {}
    for f in failures:
        for bname in BODIES.values():
            if bname in f:
                cats[bname] = cats.get(bname, 0) + 1
    for cat, count in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}")
else:
    print("\nAll tests passed!")
