"""Round 137: Lilith Sidereal Mode Deep Sweep.

Tests Mean Lilith and Osculating Lilith in sidereal mode across all
ayanamsha modes, multiple epochs, and flag combinations.
"""

from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_SIDEREAL,
    FLG_J2000,
    FLG_NONUT,
    FLG_EQUATORIAL,
    MEAN_APOG,
    OSCU_APOG,
    TRUE_NODE,
    MEAN_NODE,
)

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    MEAN_APOG: "MeanLilith",
    OSCU_APOG: "OscuLilith",
    MEAN_NODE: "MeanNode",
    TRUE_NODE: "TrueNode",
}

TEST_JDS = [
    2433282.5,  # 1950
    2440587.5,  # 1970
    2444239.5,  # 1980
    2447892.5,  # 1990
    2451545.0,  # 2000 J2000
    2455197.5,  # 2010
    2458849.5,  # 2020
    2460676.5,  # 2025
    2462502.5,  # 2030
    2466154.5,  # 2040
]

# Sidereal modes to test (most common ones)
SIDEREAL_MODES = [
    (0, "FAGAN_BRADLEY"),
    (1, "LAHIRI"),
    (2, "DELUCE"),
    (3, "RAMAN"),
    (4, "USHASHASHI"),
    (5, "KRISHNAMURTI"),
    (6, "DJWHAL_KHUL"),
    (7, "YUKTESHWAR"),
    (8, "JN_BHASIN"),
    (9, "BABYL_KUGLER1"),
    (15, "ALDEBARAN_15TAU"),
    (21, "SASSANIAN"),
    (22, "GALCENT_0SAG"),
    (27, "TRUE_CITRA"),
    (30, "SURYASIDDHANTA"),
    (40, "TRUE_REVATI"),
]

FLAG_COMBOS = {
    "SID": FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL,
    "SID+NONUT": FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL | FLG_NONUT,
    "SID+EQUAT": FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL,
}

passed = 0
failed = 0
errors = 0
total = 0
worst_cases = []

for sid_mode, sid_name in SIDEREAL_MODES:
    swe.set_sid_mode(sid_mode)
    ephem.set_sid_mode(sid_mode)

    for jd in TEST_JDS:
        for body, bname in BODIES.items():
            for flag_name, flags in FLAG_COMBOS.items():
                try:
                    se_result = swe.calc_ut(jd, body, flags)
                    le_result = ephem.calc_ut(jd, body, flags)

                    for i, comp in enumerate(
                        ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]
                    ):
                        total += 1
                        se_val = se_result[0][i]
                        le_val = le_result[0][i]
                        diff = abs(le_val - se_val)

                        if comp == "dist" or comp == "dist_spd":
                            tol = 0.001  # AU
                            is_pass = diff < tol
                        else:
                            diff_arcsec = diff * 3600.0
                            if comp in ("lon", "lat"):
                                tol = 60.0  # 1 arcmin for sidereal (known ~14" offset)
                            elif body == OSCU_APOG and "spd" in comp:
                                tol = 200.0  # OscuLilith speed known diff
                            else:
                                tol = 60.0
                            is_pass = diff_arcsec < tol

                        if is_pass:
                            passed += 1
                        else:
                            diff_arcsec = diff * 3600.0 if "dist" not in comp else diff
                            info = f"FAIL {bname:12s} {comp:8s} {sid_name:20s} JD={jd} {flag_name:12s} SE={se_val:+14.6f} LE={le_val:+14.6f} diff={diff_arcsec:.3f}"
                            print(info)
                            worst_cases.append((diff_arcsec, info))
                            failed += 1

                except Exception as e:
                    errors += 1
                    print(
                        f"ERR  {bname:12s} {sid_name:20s} JD={jd} {flag_name:12s}: {e}"
                    )

# Reset sidereal mode
swe.set_sid_mode(0)
ephem.set_sid_mode(0)

print(f"\n{'=' * 60}")
print(f"Round 137: Lilith Sidereal Mode Deep Sweep")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")

if worst_cases:
    worst_cases.sort(key=lambda x: -x[0])
    print(f"\nTop 10 worst:")
    for _, info in worst_cases[:10]:
        print(f"  {info}")
