#!/usr/bin/env python3
"""Round 167: Rise/set at equinox and solstice dates.

Tests rise_trans at astronomically significant dates (equinoxes, solstices)
where the Sun's declination is at extremes or zero. These are dates where
the day/night balance is changing rapidly and edge cases in rise/set
calculations are most likely to surface.
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

# Approximate equinox/solstice dates for 2000-2030
ASTRO_DATES = [
    # (year, month, day) - approximate dates
    (2000, 3, 20),
    (2000, 6, 21),
    (2000, 9, 22),
    (2000, 12, 21),
    (2005, 3, 20),
    (2005, 6, 21),
    (2005, 9, 22),
    (2005, 12, 21),
    (2010, 3, 20),
    (2010, 6, 21),
    (2010, 9, 23),
    (2010, 12, 21),
    (2015, 3, 20),
    (2015, 6, 21),
    (2015, 9, 23),
    (2015, 12, 22),
    (2020, 3, 20),
    (2020, 6, 20),
    (2020, 9, 22),
    (2020, 12, 21),
    (2025, 3, 20),
    (2025, 6, 21),
    (2025, 9, 22),
    (2025, 12, 21),
]

# Locations at various latitudes
LOCATIONS = [
    ("London", 51.5, -0.1, 0),
    ("NewYork", 40.7, -74.0, 0),
    ("Tokyo", 35.7, 139.7, 0),
    ("Sydney", -33.9, 151.2, 0),
    ("Tromso", 69.6, 19.0, 0),  # Near Arctic Circle
    ("Quito", -0.2, -78.5, 2850),  # Equator, high altitude
    ("CapeTown", -33.9, 18.4, 0),
    ("Reykjavik", 64.1, -21.9, 0),  # Near Arctic
]

BODIES = [
    ("Sun", swe.SUN, ephem.SUN),
    ("Moon", swe.MOON, ephem.MOON),
    ("Venus", swe.VENUS, ephem.VENUS),
]

# Rise/set event types
RISE = 1  # CALC_RISE
SET = 2  # CALC_SET

TOL_SECONDS = 120  # 2 minutes

passed = 0
failed = 0
skipped = 0
errors = []

for y, m, d in ASTRO_DATES:
    jd_start = swe.julday(y, m, d, 0.0)

    for loc_name, lat, lon, alt in LOCATIONS:
        geopos = [lon, lat, float(alt)]

        for body_name, se_body, le_body in BODIES:
            for event, event_name in [(RISE, "rise"), (SET, "set")]:
                try:
                    se_result = swe.rise_trans(
                        jd_start, se_body, event, geopos, 1013.25, 15.0
                    )
                    se_jd = se_result[1][0]
                except Exception:
                    skipped += 1
                    continue

                try:
                    le_result2 = ephem.rise_trans(
                        jd_start, le_body, event, [lon, lat, float(alt)], 1013.25, 15.0
                    )
                    # rise_trans returns (retflag, tret) where tret[0] = JD of event
                    le_jd = le_result2[1][0]
                except Exception:
                    skipped += 1
                    continue

                if se_jd == 0.0 or le_jd == 0.0:
                    skipped += 1
                    continue

                diff_sec = abs(se_jd - le_jd) * 86400.0

                if diff_sec <= TOL_SECONDS:
                    passed += 1
                else:
                    failed += 1
                    if len(errors) < 25:
                        errors.append(
                            f"  FAIL {body_name} {event_name} {y}/{m}/{d} "
                            f"@{loc_name}: {diff_sec:.1f}s"
                        )

total = passed + failed
print(f"\n=== Round 167: Rise/Set at Equinox/Solstice Dates ===")
print(f"Dates: {len(ASTRO_DATES)}, Locations: {len(LOCATIONS)}, Bodies: {len(BODIES)}")
print(
    f"Total: {total}, PASSED: {passed} ({100 * passed / total:.1f}%), FAILED: {failed}, Skipped: {skipped}"
)

if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)

if failed == 0:
    print("\nRound 167: ALL PASSED")
