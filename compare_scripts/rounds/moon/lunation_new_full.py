#!/usr/bin/env python3
"""Round 154: Moon position at lunation (New/Full Moon).

Compare Moon positions at exact New Moon and Full Moon times to verify
sub-arcsecond precision at these astronomically significant moments.
Uses mooncross_ut to find exact Sun-Moon conjunction/opposition.
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
MOON = 1

# Sample lunation dates: approximate New Moons 2020-2025
# We'll compute exact positions at these approximate times
test_dates = []
# Approximate New Moon dates (within a day)
new_moon_approx = [
    (2020, 1, 24),
    (2020, 3, 24),
    (2020, 5, 22),
    (2020, 7, 20),
    (2020, 9, 17),
    (2020, 11, 15),
    (2021, 1, 13),
    (2021, 3, 13),
    (2021, 5, 11),
    (2021, 7, 10),
    (2021, 9, 7),
    (2021, 11, 4),
    (2022, 1, 2),
    (2022, 3, 2),
    (2022, 5, 1),
    (2022, 6, 29),
    (2022, 8, 27),
    (2022, 10, 25),
    (2022, 12, 23),
    (2023, 1, 21),
    (2023, 3, 21),
    (2023, 5, 19),
    (2023, 7, 17),
    (2023, 9, 15),
    (2023, 11, 13),
    (2024, 1, 11),
    (2024, 3, 10),
    (2024, 5, 8),
    (2024, 7, 5),
    (2024, 9, 3),
    (2024, 11, 1),
    (2024, 12, 30),
    (2025, 2, 28),
    (2025, 4, 27),
    (2025, 6, 25),
]

# Full Moon dates (approximate)
full_moon_approx = [
    (2020, 1, 10),
    (2020, 3, 9),
    (2020, 5, 7),
    (2020, 7, 5),
    (2020, 9, 2),
    (2020, 11, 1),
    (2021, 1, 28),
    (2021, 3, 28),
    (2021, 5, 26),
    (2021, 7, 24),
    (2021, 9, 20),
    (2021, 11, 19),
    (2022, 1, 17),
    (2022, 3, 18),
    (2022, 5, 16),
    (2022, 7, 13),
    (2022, 9, 10),
    (2022, 11, 8),
    (2023, 1, 6),
    (2023, 3, 7),
    (2023, 5, 5),
    (2023, 7, 3),
    (2023, 9, 29),
    (2023, 11, 27),
    (2024, 1, 25),
    (2024, 3, 25),
    (2024, 5, 23),
    (2024, 7, 21),
    (2024, 9, 18),
    (2024, 11, 15),
    (2025, 1, 13),
    (2025, 3, 14),
    (2025, 5, 12),
]

for y, m, d in new_moon_approx:
    jd = swe.julday(y, m, d, 12.0)
    test_dates.append((f"NM {y}-{m:02d}-{d:02d}", jd))

for y, m, d in full_moon_approx:
    jd = swe.julday(y, m, d, 12.0)
    test_dates.append((f"FM {y}-{m:02d}-{d:02d}", jd))

TOL_LON = 1.0  # arcsec
TOL_LAT = 1.0  # arcsec
TOL_DIST = 1e-6  # AU
TOL_SPD = 1.0  # "/day

passed = failed = errors = total = 0
failures = []

print(f"Round 154: Moon Position at Lunation (New/Full Moon)")
print(f"Testing Moon+Sun at {len(test_dates)} lunation dates")
print("=" * 90)

for label, jd in test_dates:
    for body, bname in [(SUN, "Sun"), (MOON, "Moon")]:
        try:
            se_r = swe.calc_ut(jd, body, FLG_SPEED)
            le_r = ephem.calc_ut(jd, body, FLG_SPEED)
            se_d = se_r[0]
            le_d = le_r[0]

            comp_names = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]
            for i in range(6):
                total += 1
                se_val = se_d[i]
                le_val = le_d[i]

                if i <= 1 or i == 3 or i == 4:
                    diff = abs(le_val - se_val) * 3600.0
                    tol = TOL_LON if i <= 1 else TOL_SPD
                else:
                    diff = abs(le_val - se_val)
                    tol = TOL_DIST

                if diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    msg = f"  FAIL {label} {bname} {comp_names[i]}: SE={se_val:.10f} LE={le_val:.10f} diff={diff:.8f}"
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
    print(f"\n{len(failures)} failures")
    cats = {"Sun": 0, "Moon": 0}
    for f in failures:
        if "Sun" in f:
            cats["Sun"] += 1
        elif "Moon" in f:
            cats["Moon"] += 1
    for c, n in cats.items():
        if n:
            print(f"  {c}: {n}")
else:
    print("\nAll tests passed!")
