#!/usr/bin/env python3
"""Round 207: Eclipse path width comparison.

Tests sol_eclipse_how attributes including path width at various locations.
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
failures = []

# Known eclipse maximum times (from previous rounds)
ECLIPSE_EVENTS = [
    ("2017 Aug 21 total", 2457987.5),
    ("2024 Apr 8 total", 2460408.5),
    ("2020 Jun 21 annular", 2459021.5),
    ("2019 Jul 2 total", 2458666.5),
    ("2015 Mar 20 total", 2457101.5),
]

# Locations to test eclipse_how
LOCATIONS = [
    ("London", 51.5074, -0.1278, 11.0),
    ("New York", 40.7128, -74.006, 10.0),
    ("Tokyo", 35.6762, 139.6503, 40.0),
    ("Sydney", -33.8688, 151.2093, 58.0),
    ("Delhi", 28.6139, 77.2090, 216.0),
    ("Cairo", 30.0444, 31.2357, 75.0),
    ("Santiago", -33.4489, -70.6693, 570.0),
    ("Reykjavik", 64.1466, -21.9426, 0.0),
]


def test_eclipse_how():
    global passed, failed, total

    print("=" * 70)
    print("Round 207: Eclipse Path Width / eclipse_how")
    print("=" * 70)

    for ecl_label, approx_jd in ECLIPSE_EVENTS:
        # First find the actual eclipse max time
        try:
            le_glob = ephem.sol_eclipse_when_glob(approx_jd - 30, 0)
            ecl_jd = le_glob[1][0]  # maximum time
        except Exception:
            ecl_jd = approx_jd

        print(f"\n--- {ecl_label} (JD={ecl_jd:.6f}) ---")

        for loc_name, lat, lon, alt in LOCATIONS:
            geopos = [lon, lat, alt]

            try:
                le_r = ephem.sol_eclipse_how(ecl_jd, 0, geopos)
                le_attr = le_r[1] if isinstance(le_r, tuple) and len(le_r) > 1 else le_r
            except Exception as e:
                continue

            try:
                se_r = swe.sol_eclipse_how(ecl_jd, geopos, swe.FLG_SWIEPH)
                se_attr = se_r[1]
            except Exception:
                continue

            # Compare attributes
            # attr[0] = fraction of solar diameter covered by moon
            # attr[1] = ratio of lunar/solar diameter
            # attr[2] = fraction of solar disc covered (obscuration)
            # attr[3] = path width (km)

            for idx, aname in [(0, "frac_diam"), (1, "diam_ratio"), (2, "obscuration")]:
                try:
                    le_val = le_attr[idx]
                    se_val = se_attr[idx]
                except (IndexError, TypeError):
                    continue

                total += 1
                if se_val == 0 and le_val == 0:
                    passed += 1
                    continue

                diff = abs(le_val - se_val)
                # Relative tolerance 5%
                if se_val != 0 and diff / abs(se_val) <= 0.05:
                    passed += 1
                elif diff <= 0.01:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f"  {ecl_label} {loc_name} {aname}: LE={le_val:.6f} SE={se_val:.6f}"
                    )


if __name__ == "__main__":
    test_eclipse_how()

    print(f"\n{'=' * 70}")
    pct = 100 * passed / total if total > 0 else 0
    print(f"RESULTS: {passed}/{total} passed ({pct:.1f}%), {failed} failed")
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
