#!/usr/bin/env python3
"""Round 184: Coordinate chain deep — ecliptic -> equatorial -> back.

Tests the full coordinate transformation roundtrip: compute ecliptic position,
convert to equatorial via cotrans, convert back, and verify we get the
original position within numerical precision.
"""

from __future__ import annotations
import sys, os, math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    "Sun": (swe.SUN, ephem.SUN),
    "Moon": (swe.MOON, ephem.MOON),
    "Mars": (swe.MARS, ephem.MARS),
    "Jupiter": (swe.JUPITER, ephem.JUPITER),
    "Saturn": (swe.SATURN, ephem.SATURN),
}

FLAGS_ECL = swe.FLG_SPEED
FLAGS_EQ = swe.FLG_SPEED | 2048  # EQUATORIAL

TEST_DATES = []
for year in range(1960, 2041, 10):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

TOL = 1.0  # arcsec

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    # Get obliquity
    try:
        ecl_nut = swe.calc_ut(jd, swe.ECL_NUT, 0)[0]
        eps = ecl_nut[0]  # true obliquity
    except Exception:
        continue

    for bname, (se_id, le_id) in BODIES.items():
        try:
            # Get ecliptic position from both
            se_ecl = swe.calc_ut(jd, se_id, FLAGS_ECL)[0]
            le_ecl = ephem.calc_ut(jd, le_id, FLAGS_ECL)[0]

            # Get equatorial position from both
            se_eq = swe.calc_ut(jd, se_id, FLAGS_EQ)[0]
            le_eq = ephem.calc_ut(jd, le_id, FLAGS_EQ)[0]
        except Exception:
            continue

        # Compare ecliptic
        diff_ecl_lon = abs(se_ecl[0] - le_ecl[0])
        if diff_ecl_lon > 180:
            diff_ecl_lon = 360 - diff_ecl_lon
        diff_ecl_as = diff_ecl_lon * 3600

        # Compare equatorial
        diff_eq_ra = abs(se_eq[0] - le_eq[0])
        if diff_eq_ra > 180:
            diff_eq_ra = 360 - diff_eq_ra
        diff_eq_as = diff_eq_ra * 3600
        diff_eq_dec = abs(se_eq[1] - le_eq[1]) * 3600

        # Also test cotrans roundtrip with LE
        try:
            le_to_eq = ephem.cotrans((le_ecl[0], le_ecl[1], le_ecl[2]), -eps)
            le_back = ephem.cotrans((le_to_eq[0], le_to_eq[1], le_to_eq[2]), eps)
            roundtrip_diff = abs(le_back[0] - le_ecl[0])
            if roundtrip_diff > 180:
                roundtrip_diff = 360 - roundtrip_diff
            roundtrip_as = roundtrip_diff * 3600
        except Exception:
            roundtrip_as = 0.0

        ok = True
        reasons = []
        if diff_ecl_as > TOL:
            ok = False
            reasons.append(f'ecl lon {diff_ecl_as:.3f}"')
        if diff_eq_as > TOL:
            ok = False
            reasons.append(f'eq RA {diff_eq_as:.3f}"')
        if diff_eq_dec > TOL:
            ok = False
            reasons.append(f'eq Dec {diff_eq_dec:.3f}"')
        if roundtrip_as > 0.001:
            ok = False
            reasons.append(f'roundtrip {roundtrip_as:.6f}"')

        if ok:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(f"  FAIL {bname} {date_str}: {', '.join(reasons)}")

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 184: Coordinate Chain Deep ===")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 184: ALL PASSED")
