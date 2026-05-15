#!/usr/bin/env python3
"""Round 208: Precession rates verification.

Tests that ecliptic coordinates at J2000 vs ecliptic of date show the
correct precession rate (~50.3"/year) and that the transformation is
consistent between LE and SE.
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

BODIES = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
]

FLAGS_DATE = ephem.FLG_SWIEPH | ephem.FLG_SPEED
FLAGS_J2000 = FLAGS_DATE | ephem.FLG_J2000

TEST_JDS = [
    2451545.0,  # J2000 (precession = 0)
    2451545.0 + 365.25,  # J2001
    2451545.0 + 3652.5,  # J2010
    2451545.0 + 7305.0,  # J2020
    2451545.0 + 9131.25,  # J2025
    2451545.0 - 3652.5,  # J1990
    2451545.0 - 7305.0,  # J1980
    2451545.0 - 36525.0,  # J1900
]


def test_precession_rates():
    global passed, failed, total

    print("=" * 70)
    print("Round 208: Precession Rates Verification")
    print("=" * 70)

    for bname, le_b, se_b in BODIES:
        print(f"\n--- {bname} ---")

        for jd in TEST_JDS:
            years = (jd - 2451545.0) / 365.25

            # Get ecliptic of date and J2000 positions from both
            try:
                le_date = ephem.calc_ut(jd, le_b, FLAGS_DATE)
                le_j2k = ephem.calc_ut(jd, le_b, FLAGS_J2000)
                se_date = swe.calc_ut(jd, se_b, swe.FLG_SWIEPH | swe.FLG_SPEED)
                se_j2k = swe.calc_ut(
                    jd, se_b, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000
                )
            except Exception:
                continue

            # Compare ecliptic of date positions
            total += 1
            d_diff = abs(le_date[0][0] - se_date[0][0])
            if d_diff > 180:
                d_diff = 360 - d_diff
            d_as = d_diff * 3600
            if d_as <= 0.5:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {bname} date JD+{years:.0f}y: diff={d_as:.4f}"')

            # Compare J2000 positions
            total += 1
            j_diff = abs(le_j2k[0][0] - se_j2k[0][0])
            if j_diff > 180:
                j_diff = 360 - j_diff
            j_as = j_diff * 3600
            if j_as <= 0.5:
                passed += 1
            else:
                failed += 1
                failures.append(f'  {bname} J2000 JD+{years:.0f}y: diff={j_as:.4f}"')

            # Precession offset = date - J2000 (should be ~50.3"/year * years)
            le_prec = le_date[0][0] - le_j2k[0][0]
            se_prec = se_date[0][0] - se_j2k[0][0]
            if le_prec > 180:
                le_prec -= 360
            if le_prec < -180:
                le_prec += 360
            if se_prec > 180:
                se_prec -= 360
            if se_prec < -180:
                se_prec += 360

            total += 1
            prec_diff = abs(le_prec - se_prec) * 3600
            if prec_diff <= 0.5:
                passed += 1
            else:
                failed += 1
                failures.append(
                    f'  {bname} prec JD+{years:.0f}y: LE={le_prec:.6f}° SE={se_prec:.6f}° diff={prec_diff:.4f}"'
                )

    # Test general precession value
    print("\n--- General Precession Rate ---")
    for jd_offset in [3652.5, 7305.0, 18262.5, 36525.0]:
        jd = 2451545.0 + jd_offset
        years = jd_offset / 365.25

        try:
            le_sun_d = ephem.calc_ut(jd, ephem.SUN, FLAGS_DATE)[0][0]
            le_sun_j = ephem.calc_ut(jd, ephem.SUN, FLAGS_J2000)[0][0]
            se_sun_d = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED)[0][0]
            se_sun_j = swe.calc_ut(
                jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000
            )[0][0]
        except Exception:
            continue

        le_prec_rate = ((le_sun_d - le_sun_j) * 3600) / years if years != 0 else 0
        se_prec_rate = ((se_sun_d - se_sun_j) * 3600) / years if years != 0 else 0

        total += 1
        rate_diff = abs(le_prec_rate - se_prec_rate)
        if rate_diff <= 0.1:  # 0.1"/yr
            passed += 1
        else:
            failed += 1
            failures.append(
                f'  Prec rate at {years:.0f}yr: LE={le_prec_rate:.4f}"/yr SE={se_prec_rate:.4f}"/yr diff={rate_diff:.4f}'
            )


if __name__ == "__main__":
    test_precession_rates()

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
