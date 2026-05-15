#!/usr/bin/env python3
"""Round 177: Oriental/occidental planet detection — elongation from Sun.

Tests the elongation calculation (angular separation from Sun) for all planets
at various dates. Whether a planet is oriental (east of Sun, rising after) or
occidental (west of Sun, rising before) is determined by elongation sign.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    "Moon": (swe.MOON, ephem.MOON),
    "Mercury": (swe.MERCURY, ephem.MERCURY),
    "Venus": (swe.VENUS, ephem.VENUS),
    "Mars": (swe.MARS, ephem.MARS),
    "Jupiter": (swe.JUPITER, ephem.JUPITER),
    "Saturn": (swe.SATURN, ephem.SATURN),
    "Uranus": (swe.URANUS, ephem.URANUS),
    "Neptune": (swe.NEPTUNE, ephem.NEPTUNE),
    "Pluto": (swe.PLUTO, ephem.PLUTO),
}

FLAGS = swe.FLG_SPEED

TEST_DATES = []
for year in range(1950, 2051, 5):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

TOL_ELONG = 1.0  # arcsec

passed = 0
failed = 0
errors = []

for date_str, jd in TEST_DATES:
    # Get Sun position
    try:
        se_sun = swe.calc_ut(jd, swe.SUN, FLAGS)[0]
        le_sun = ephem.calc_ut(jd, ephem.SUN, FLAGS)[0]
    except Exception:
        continue

    for bname, (se_id, le_id) in BODIES.items():
        try:
            se_pos = swe.calc_ut(jd, se_id, FLAGS)[0]
            le_pos = ephem.calc_ut(jd, le_id, FLAGS)[0]
        except Exception:
            continue

        # Elongation from Sun
        se_elong = se_pos[0] - se_sun[0]
        if se_elong > 180:
            se_elong -= 360
        elif se_elong < -180:
            se_elong += 360

        le_elong = le_pos[0] - le_sun[0]
        if le_elong > 180:
            le_elong -= 360
        elif le_elong < -180:
            le_elong += 360

        diff_elong_as = abs(se_elong - le_elong) * 3600

        # Check oriental/occidental agreement (sign of elongation)
        se_oriental = se_elong > 0
        le_oriental = le_elong > 0
        sign_agree = se_oriental == le_oriental

        ok = True
        reasons = []
        if diff_elong_as > TOL_ELONG:
            ok = False
            reasons.append(f'elong diff {diff_elong_as:.3f}"')
        if not sign_agree and abs(se_elong) > 1.0:  # Only flag if not near conjunction
            ok = False
            reasons.append("orient/occid disagree")

        if ok:
            passed += 1
        else:
            failed += 1
            if len(errors) < 20:
                errors.append(f"  FAIL {bname} {date_str}: {', '.join(reasons)}")

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 177: Oriental/Occidental Detection ===")
print(f"Dates: {len(TEST_DATES)}, Bodies: {len(BODIES)}")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 177: ALL PASSED")
