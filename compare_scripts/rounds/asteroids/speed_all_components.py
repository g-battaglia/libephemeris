"""Round 136: Asteroid Speed Accuracy Deep Sweep.

Tests all 6 position/speed components for Ceres, Pallas, Juno, Vesta, Chiron
across many epochs within their valid ephemeris range, with multiple flag combos.
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
    FLG_HELCTR,
    FLG_TRUEPOS,
    FLG_J2000,
    FLG_NONUT,
    FLG_NOABERR,
    FLG_EQUATORIAL,
    FLG_SWIEPH,
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
)

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    CERES: "Ceres",
    PALLAS: "Pallas",
    JUNO: "Juno",
    VESTA: "Vesta",
    CHIRON: "Chiron",
}

# Valid range: ~1930-2070 for asteroids in DE440
TEST_JDS = [
    2425977.5,  # 1930-01-01
    2429616.5,  # 1940-01-01
    2433282.5,  # 1950-01-01
    2436934.5,  # 1960-01-01
    2440587.5,  # 1970-01-01
    2444239.5,  # 1980-01-01
    2447892.5,  # 1990-01-01
    2451545.0,  # 2000-01-01.5
    2451911.5,  # 2001-01-01
    2453371.5,  # 2005-01-01
    2455197.5,  # 2010-01-01
    2456293.5,  # 2013-01-01
    2457023.5,  # 2015-01-01
    2458119.5,  # 2018-01-01
    2458849.5,  # 2020-01-01
    2459580.5,  # 2022-01-01
    2460310.5,  # 2024-01-01
    2460676.5,  # 2025-01-01
    2461041.5,  # 2026-01-01
    2462502.5,  # 2030-01-01
    2464328.5,  # 2035-01-01
    2466154.5,  # 2040-01-01
    2467980.5,  # 2045-01-01
    2469807.5,  # 2050-01-01
    2473459.5,  # 2060-01-01
]

COMP_NAMES = ["lon", "lat", "dist", "lon_speed", "lat_speed", "dist_speed"]

FLAG_COMBOS = {
    "default": FLG_SWIEPH | FLG_SPEED,
    "J2000": FLG_SWIEPH | FLG_SPEED | FLG_J2000,
    "NONUT": FLG_SWIEPH | FLG_SPEED | FLG_NONUT,
    "TRUEPOS": FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS,
    "NOABERR": FLG_SWIEPH | FLG_SPEED | FLG_NOABERR,
    "HELCTR": FLG_SWIEPH | FLG_SPEED | FLG_HELCTR,
    "EQUATORIAL": FLG_SWIEPH | FLG_SPEED | FLG_EQUATORIAL,
}

passed = 0
failed = 0
errors = 0
total = 0
worst_cases = []

for jd in TEST_JDS:
    for body, bname in BODIES.items():
        for flag_name, flags in FLAG_COMBOS.items():
            try:
                se_result = swe.calc_ut(jd, body, flags)
                le_result = ephem.calc_ut(jd, body, flags)

                for i, comp in enumerate(COMP_NAMES):
                    total += 1
                    se_val = se_result[0][i]
                    le_val = le_result[0][i]
                    diff = abs(le_val - se_val)

                    if i < 3:  # position components
                        if comp == "dist":
                            # Distance in AU — tolerance ~10 micro-AU
                            tol = 1e-5
                            is_pass = diff < tol
                            diff_str = f"{diff:.2e} AU"
                        else:
                            # lon/lat in degrees — tolerance 5"
                            diff_arcsec = diff * 3600.0
                            tol = 5.0  # arcsec
                            is_pass = diff_arcsec < tol
                            diff_str = f'{diff_arcsec:.4f}"'
                    else:  # speed components
                        diff_arcsec = diff * 3600.0
                        if comp == "dist_speed":
                            tol = 0.01  # AU/day
                            is_pass = diff < tol
                            diff_str = f"{diff:.2e} AU/day"
                        else:
                            # Speed tolerance: 2" /day or 1% of magnitude
                            speed_mag = max(abs(se_val), abs(le_val))
                            rel_tol = speed_mag * 3600.0 * 0.01
                            abs_tol = 2.0  # "/day
                            tol = max(abs_tol, rel_tol)
                            is_pass = diff_arcsec < tol
                            diff_str = f'{diff_arcsec:.4f}"/day'

                    if is_pass:
                        passed += 1
                    else:
                        failed += 1
                        info = f"FAIL {bname:8s} {comp:10s} JD={jd} {flag_name:12s} SE={se_val:+14.8f} LE={le_val:+14.8f} diff={diff_str}"
                        print(info)
                        worst_cases.append((diff, info))

            except Exception as e:
                errors += 1
                if "Invalid Time" not in str(e):
                    print(f"ERR  {bname:8s} JD={jd} {flag_name:12s}: {e}")

print(f"\n{'=' * 60}")
print(f"Round 136: Asteroid Speed Accuracy Deep Sweep")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")

if worst_cases:
    worst_cases.sort(key=lambda x: -x[0])
    print(f"\nTop 10 worst cases:")
    for _, info in worst_cases[:10]:
        print(f"  {info}")
