#!/usr/bin/env python3
"""Round 189: Rise/set at extreme polar latitudes (75°-85°).

Tests rise/transit/set calculations at high latitudes where circumpolar
and never-rise conditions create edge cases. Tests Sun, Moon, and planets
at latitudes 75°, 78°, 80°, 82°, 85° N and S.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
skipped = 0
failures = []

# Bodies to test
BODIES = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Moon", ephem.MOON, swe.MOON),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
]

# Extreme latitudes
LATITUDES = [75.0, 78.0, 80.0, 82.0, 85.0, -75.0, -80.0, -85.0]

# Test dates across seasons
TEST_DATES = [
    ("Winter solstice 2020", 2459204.5),  # Dec 21, 2020
    ("Spring equinox 2020", 2458930.5),  # Mar 20, 2020
    ("Summer solstice 2020", 2459022.5),  # Jun 20, 2020
    ("Autumn equinox 2020", 2459115.5),  # Sep 22, 2020
    ("Winter solstice 2015", 2457377.5),  # Dec 22, 2015
    ("Summer solstice 2015", 2457194.5),  # Jun 21, 2015
]

# Rise/set event types
CALC_RISE = 1
CALC_SET = 2
CALC_MTRANSIT = 4
CALC_ITRANSIT = 8


def test_rise_set_polar():
    global passed, failed, total, skipped

    print("=" * 70)
    print("Round 189: Rise/Set at Extreme Polar Latitudes (75°-85°)")
    print("=" * 70)

    for date_label, jd in TEST_DATES:
        print(f"\n--- {date_label} (JD {jd}) ---")

        for body_name, le_body, se_body in BODIES:
            for lat in LATITUDES:
                lon = 25.0  # Fixed longitude
                alt = 0.0
                pressure = 1013.25
                temp = 10.0

                for rsmi, event_name in [
                    (CALC_RISE, "rise"),
                    (CALC_SET, "set"),
                    (CALC_MTRANSIT, "transit"),
                    (CALC_ITRANSIT, "itransit"),
                ]:
                    label = f"{body_name} {event_name} lat={lat}"

                    # pyswisseph
                    try:
                        se_result = swe.rise_trans(
                            jd, se_body, rsmi, [lon, lat, alt], pressure, temp
                        )
                        se_jd = se_result[1][0]
                        se_ok = True
                    except Exception as e:
                        se_ok = False
                        se_err = str(e)

                    # libephemeris
                    try:
                        le_result = ephem.rise_trans(
                            jd, le_body, rsmi, [lon, lat, alt], pressure, temp
                        )
                        le_jd = le_result[1][0]
                        le_ok = True
                    except Exception as e:
                        le_ok = False
                        le_err = str(e)

                    total += 1

                    if not se_ok and not le_ok:
                        # Both error — agree it's impossible
                        passed += 1
                        continue

                    if se_ok and not le_ok:
                        # SE found event but LE didn't
                        failed += 1
                        failures.append(
                            f"  {date_label} {label}: SE found JD={se_jd:.6f}, LE error: {le_err[:60]}"
                        )
                        continue

                    if not se_ok and le_ok:
                        # LE found event but SE didn't
                        failed += 1
                        failures.append(
                            f"  {date_label} {label}: LE found JD={le_jd:.6f}, SE error: {se_err[:60]}"
                        )
                        continue

                    # Both found — compare
                    if se_jd == 0.0 and le_jd == 0.0:
                        passed += 1
                        continue

                    if se_jd == 0.0 or le_jd == 0.0:
                        failed += 1
                        failures.append(
                            f"  {date_label} {label}: SE_JD={se_jd:.6f} LE_JD={le_jd:.6f} (one is zero)"
                        )
                        continue

                    diff_min = abs(le_jd - se_jd) * 1440  # minutes

                    # Tolerance: 5 minutes for rise/set at polar, 2 min for transit
                    tol = 2.0 if "transit" in event_name else 5.0
                    if diff_min <= tol:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f"  {date_label} {label}: diff={diff_min:.2f} min (tol {tol})"
                        )


if __name__ == "__main__":
    test_rise_set_polar()

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
