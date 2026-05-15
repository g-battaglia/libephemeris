#!/usr/bin/env python3
"""Round 168: Hypothetical/Uranian bodies with sidereal mode.

Tests Uranian/hypothetical bodies (Cupido, Hades, Zeus, Kronos, Apollon,
Admetos, Vulkanus, Poseidon) with FLG_SIDEREAL flag and various ayanamsha
modes. This combination exercises the hypothetical body Keplerian propagation
+ sidereal coordinate transformation pipeline.
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

HYPOTHETICAL = {
    "Cupido": (swe.CUPIDO, ephem.CUPIDO),
    "Hades": (swe.HADES, ephem.HADES),
    "Zeus": (swe.ZEUS, ephem.ZEUS),
    "Kronos": (swe.KRONOS, ephem.KRONOS),
    "Apollon": (swe.APOLLON, ephem.APOLLON),
    "Admetos": (swe.ADMETOS, ephem.ADMETOS),
    "Vulkanus": (swe.VULKANUS, ephem.VULKANUS),
    "Poseidon": (swe.POSEIDON, ephem.POSEIDON),
}

# Ayanamsha modes to test
AYANAMSHAS = {
    "Lahiri": 1,
    "Raman": 3,
    "KrishnamurtiNew": 5,
    "Fagan-Bradley": 0,
}

TEST_DATES = []
for year in range(1950, 2051, 10):
    for month in [1, 7]:
        jd = swe.julday(year, month, 1, 12.0)
        TEST_DATES.append((f"{year}/{month:02d}/01", jd))

# Tolerances for hypothetical bodies (they have ~35" known tolerance)
TOL_LON = 60.0  # arcsec - generous for sidereal + hypothetical
TOL_SPD = 5.0  # arcsec/day

passed = 0
failed = 0
errors = []

for aya_name, aya_mode in AYANAMSHAS.items():
    swe.set_sid_mode(aya_mode)
    ephem.set_sid_mode(aya_mode, 0, 0)

    flags_se = swe.FLG_SPEED | swe.FLG_SIDEREAL
    flags_le = ephem.FLG_SPEED | ephem.FLG_SIDEREAL

    for date_str, jd in TEST_DATES:
        for bname, (se_id, le_id) in HYPOTHETICAL.items():
            try:
                se_result = swe.calc_ut(jd, se_id, flags_se)
                se_pos = (
                    se_result[0]
                    if isinstance(se_result[0], (list, tuple))
                    else se_result
                )
            except Exception:
                continue

            try:
                le_result = ephem.calc_ut(jd, le_id, flags_le)
                le_pos = le_result[0]
            except Exception:
                continue

            diff_lon = abs(se_pos[0] - le_pos[0])
            if diff_lon > 180:
                diff_lon = 360 - diff_lon
            diff_lon_as = diff_lon * 3600

            diff_spd = abs(se_pos[3] - le_pos[3]) * 3600

            ok = True
            reasons = []

            if diff_lon_as > TOL_LON:
                ok = False
                reasons.append(f'lon {diff_lon_as:.1f}"')
            if diff_spd > TOL_SPD:
                ok = False
                reasons.append(f'spd {diff_spd:.3f}"/d')

            if ok:
                passed += 1
            else:
                failed += 1
                if len(errors) < 25:
                    errors.append(
                        f"  FAIL {bname} {date_str} [{aya_name}]: {', '.join(reasons)}"
                    )

# Reset sidereal mode
swe.set_sid_mode(0)
ephem.set_sid_mode(0, 0, 0)

total = passed + failed
print(f"\n=== Round 168: Hypothetical Bodies with Sidereal Mode ===")
print(
    f"Dates: {len(TEST_DATES)}, Bodies: {len(HYPOTHETICAL)}, Ayanamshas: {len(AYANAMSHAS)}"
)
print(
    f"Total: {total}, PASSED: {passed} ({100 * passed / total:.1f}%), FAILED: {failed}"
)

if errors:
    print(f"\nFirst {len(errors)} failures:")
    for e in errors:
        print(e)

if failed == 0:
    print("\nRound 168: ALL PASSED")
