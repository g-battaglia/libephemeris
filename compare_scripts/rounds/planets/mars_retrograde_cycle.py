"""Round 141: Mars Retrograde Cycle Deep.

Tests Mars positions at fine intervals during retrograde periods,
verifying longitude reversal, speed sign change, and station timing.
"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import FLG_SPEED, FLG_SWIEPH, MARS

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

flags = FLG_SWIEPH | FLG_SPEED

# Mars retrograde periods (approximate JD ranges around station/retrograde)
# Each is (start_jd, end_jd, label) covering ~4 months around opposition
RETRO_PERIODS = [
    (2451800.0, 2451950.0, "2000-2001"),  # Jun-Nov 2000
    (2452920.0, 2453070.0, "2003"),  # Jul-Dec 2003
    (2454040.0, 2454190.0, "2005-2006"),  # Oct 2005 - Mar 2006
    (2455460.0, 2455610.0, "2010"),  # Dec 2009 - May 2010
    (2456330.0, 2456480.0, "2012"),  # Jan-Apr 2012
    (2457100.0, 2457250.0, "2014"),  # Mar-Jun 2014
    (2457660.0, 2457810.0, "2016"),  # Apr-Jul 2016
    (2458400.0, 2458550.0, "2018"),  # Jul-Sep 2018
    (2458750.0, 2458900.0, "2019-2020"),  # Sep 2019 - Jan 2020
    (2459440.0, 2459590.0, "2022"),  # Oct 2022 - Jan 2023
    (2460250.0, 2460400.0, "2024"),  # Dec 2024 - Feb 2025
]

passed = 0
failed = 0
total = 0

for start_jd, end_jd, label in RETRO_PERIODS:
    # Sample every 2 days across the period
    jd = start_jd
    se_stations = []
    le_stations = []
    prev_se_spd = None
    prev_le_spd = None

    while jd <= end_jd:
        total += 1
        try:
            se_result = swe.calc_ut(jd, MARS, flags)
            le_result = ephem.calc_ut(jd, MARS, flags)

            se_lon, se_lat = se_result[0][0], se_result[0][1]
            le_lon, le_lat = le_result[0][0], le_result[0][1]
            se_spd = se_result[0][3]
            le_spd = le_result[0][3]

            lon_diff = abs(le_lon - se_lon) * 3600.0
            lat_diff = abs(le_lat - se_lat) * 3600.0
            spd_diff = abs(le_spd - se_spd) * 3600.0

            # Track stations (speed sign changes)
            if prev_se_spd is not None and prev_se_spd * se_spd < 0:
                se_stations.append(jd)
            if prev_le_spd is not None and prev_le_spd * le_spd < 0:
                le_stations.append(jd)

            prev_se_spd = se_spd
            prev_le_spd = le_spd

            # Tolerances
            lon_tol = 1.0  # 1" position
            spd_tol = 2.0  # 2"/day speed

            if lon_diff < lon_tol and spd_diff < spd_tol:
                passed += 1
            else:
                failed += 1
                if lon_diff >= lon_tol:
                    print(
                        f'FAIL Mars lon {label} JD={jd:.1f} SE={se_lon:.8f} LE={le_lon:.8f} diff={lon_diff:.4f}"'
                    )
                if spd_diff >= spd_tol:
                    print(
                        f'FAIL Mars spd {label} JD={jd:.1f} SE={se_spd:.8f} LE={le_spd:.8f} diff={spd_diff:.4f}"/day'
                    )

        except Exception as e:
            failed += 1
            print(f"ERR  Mars {label} JD={jd:.1f}: {e}")

        jd += 2.0

    # Check station count matches
    total += 1
    if len(se_stations) == len(le_stations):
        passed += 1
    else:
        failed += 1
        print(
            f"FAIL Mars station count {label}: SE={len(se_stations)} LE={len(le_stations)}"
        )

print(f"\n{'=' * 60}")
print(f"Round 141: Mars Retrograde Cycle Deep")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
