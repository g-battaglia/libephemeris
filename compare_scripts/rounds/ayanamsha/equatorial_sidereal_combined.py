"""Round 139: All Planets Equatorial+Sidereal Combined.

Tests the combination of FLG_EQUATORIAL + FLG_SIDEREAL for all planets.
This is a tricky flag combination where ayanamsha must be applied correctly
to RA/Dec coordinates.
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
    FLG_EQUATORIAL,
    FLG_NONUT,
    FLG_J2000,
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
)

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    SUN: "Sun",
    MOON: "Moon",
    MERCURY: "Mercury",
    VENUS: "Venus",
    MARS: "Mars",
    JUPITER: "Jupiter",
    SATURN: "Saturn",
    URANUS: "Uranus",
    NEPTUNE: "Neptune",
    PLUTO: "Pluto",
    CHIRON: "Chiron",
    CERES: "Ceres",
    PALLAS: "Pallas",
    JUNO: "Juno",
    VESTA: "Vesta",
    MEAN_NODE: "MeanNode",
    TRUE_NODE: "TrueNode",
    MEAN_APOG: "MeanApog",
    OSCU_APOG: "OscuApog",
}

TEST_JDS = [
    2433282.5,
    2440587.5,
    2444239.5,
    2447892.5,
    2451545.0,
    2455197.5,
    2458849.5,
    2460676.5,
    2462502.5,
    2466154.5,
]

SIDEREAL_MODES = [0, 1, 3, 5, 27]  # Fagan, Lahiri, Raman, Krishnamurti, TrueCitra

FLAG_COMBOS = {
    "SID+EQ": FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL | FLG_EQUATORIAL,
    "SID+EQ+NONUT": FLG_SWIEPH
    | FLG_SPEED
    | FLG_SIDEREAL
    | FLG_EQUATORIAL
    | FLG_NONUT,
}

COMP_NAMES = ["RA", "Dec", "dist", "RA_spd", "Dec_spd", "dist_spd"]

passed = 0
failed = 0
errors = 0
total = 0

for sid_mode in SIDEREAL_MODES:
    swe.set_sid_mode(sid_mode)
    ephem.set_sid_mode(sid_mode)

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

                        if comp == "dist" or comp == "dist_spd":
                            tol = 0.001
                            is_pass = diff < tol
                        else:
                            diff_arcsec = diff * 3600.0
                            if body in (
                                MEAN_NODE,
                                TRUE_NODE,
                                MEAN_APOG,
                                OSCU_APOG,
                            ):
                                tol = 120.0  # 2 arcmin for analytical bodies
                            elif "spd" in comp and body == OSCU_APOG:
                                tol = 200.0
                            else:
                                tol = 60.0  # 1 arcmin (includes ~14" sidereal offset)
                            is_pass = diff_arcsec < tol

                        if is_pass:
                            passed += 1
                        else:
                            diff_arcsec = diff * 3600.0
                            info = f'FAIL {bname:10s} {comp:8s} sid={sid_mode} JD={jd} {flag_name} SE={se_val:+14.6f} LE={le_val:+14.6f} diff={diff_arcsec:.3f}"'
                            print(info)
                            failed += 1

                except Exception as e:
                    errors += 1
                    if "Invalid Time" not in str(e):
                        print(
                            f"ERR  {bname:10s} sid={sid_mode} JD={jd} {flag_name}: {e}"
                        )

swe.set_sid_mode(0)
ephem.set_sid_mode(0)

print(f"\n{'=' * 60}")
print(f"Round 139: All Planets Equatorial+Sidereal Combined")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")
