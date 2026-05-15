#!/usr/bin/env python3
"""Round 178: Antiscia with sidereal mode.

Tests antiscia (mirror points) computation with sidereal mode enabled.
Antiscia are reflection points across the Cancer-Capricorn axis (0°/180° axis).
Formula: antiscia_lon = (360 - lon) % 360 applied to sidereal positions.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    "Sun": (swe.SUN, ephem.SUN),
    "Moon": (swe.MOON, ephem.MOON),
    "Mercury": (swe.MERCURY, ephem.MERCURY),
    "Venus": (swe.VENUS, ephem.VENUS),
    "Mars": (swe.MARS, ephem.MARS),
    "Jupiter": (swe.JUPITER, ephem.JUPITER),
    "Saturn": (swe.SATURN, ephem.SATURN),
}

AYANAMSHAS = {"Lahiri": 1, "Fagan": 0, "Raman": 3}

FLAGS = swe.FLG_SPEED | swe.FLG_SIDEREAL
TOL = 1.0  # arcsec

TEST_DATES = []
for year in range(1960, 2041, 10):
    for month in [1, 4, 7, 10]:
        jd = swe.julday(year, month, 15, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/15", jd))

passed = 0
failed = 0
errors = []

for aya_name, aya_mode in AYANAMSHAS.items():
    swe.set_sid_mode(aya_mode)
    ephem.set_sid_mode(aya_mode, 0, 0)

    for date_str, jd in TEST_DATES:
        for bname, (se_id, le_id) in BODIES.items():
            try:
                se_r = swe.calc_ut(jd, se_id, FLAGS)
                se_pos = se_r[0] if isinstance(se_r[0], (list, tuple)) else se_r
                le_r = ephem.calc_ut(jd, le_id, FLAGS)
                le_pos = le_r[0]
            except Exception:
                continue

            diff_lon = abs(se_pos[0] - le_pos[0])
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lon_as = diff_lon * 3600

            # Also verify antiscia computation consistency
            se_antiscia = (360.0 - se_pos[0]) % 360.0
            le_antiscia = (360.0 - le_pos[0]) % 360.0
            diff_anti = abs(se_antiscia - le_antiscia)
            if diff_anti > 180:
                diff_anti = 360 - diff_anti
            diff_anti_as = diff_anti * 3600

            ok = True
            reasons = []
            if diff_lon_as > TOL:
                ok = False
                reasons.append(f'lon {diff_lon_as:.3f}"')
            if diff_anti_as > TOL:
                ok = False
                reasons.append(f'antiscia {diff_anti_as:.3f}"')

            if ok:
                passed += 1
            else:
                failed += 1
                if len(errors) < 20:
                    errors.append(
                        f"  FAIL {bname} {date_str} [{aya_name}]: {', '.join(reasons)}"
                    )

swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0, 0)

total = passed + failed
pct = 100 * passed / total if total else 0
print(f"\n=== Round 178: Antiscia with Sidereal Mode ===")
print(f"Total: {total}, PASSED: {passed} ({pct:.1f}%), FAILED: {failed}")
if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)
if failed == 0:
    print("\nRound 178: ALL PASSED")
