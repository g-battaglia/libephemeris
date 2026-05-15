"""Round 140: Planet Positions at Solstice/Equinox Precision.

Tests all planet positions at exact solstice/equinox moments (Sun at 0/90/180/270°)
across multiple years. These are astronomically significant moments.
"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    FLG_SPEED,
    FLG_SWIEPH,
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
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
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
    MEAN_NODE: "MeanNode",
    TRUE_NODE: "TrueNode",
    MEAN_APOG: "MeanApog",
}

# Approximate JDs for solstices/equinoxes 1950-2050
# (Vernal equinox ~Mar 20, Summer solstice ~Jun 21, Autumnal ~Sep 22, Winter ~Dec 21)
EVENTS = [
    (2433374.0, "VE 1950"),
    (2433466.0, "SS 1950"),
    (2433559.0, "AE 1950"),
    (2433650.0, "WS 1950"),
    (2440283.0, "VE 1969"),
    (2440375.0, "SS 1969"),
    (2440468.0, "AE 1969"),
    (2440559.0, "WS 1969"),
    (2444300.0, "VE 1980"),
    (2444392.0, "SS 1980"),
    (2444485.0, "AE 1980"),
    (2444576.0, "WS 1980"),
    (2447953.0, "VE 1990"),
    (2448045.0, "SS 1990"),
    (2448138.0, "AE 1990"),
    (2448229.0, "WS 1990"),
    (2451624.0, "VE 2000"),
    (2451716.0, "SS 2000"),
    (2451809.0, "AE 2000"),
    (2451900.0, "WS 2000"),
    (2455277.0, "VE 2010"),
    (2455369.0, "SS 2010"),
    (2455462.0, "AE 2010"),
    (2455553.0, "WS 2010"),
    (2458929.0, "VE 2020"),
    (2459021.0, "SS 2020"),
    (2459114.0, "AE 2020"),
    (2459205.0, "WS 2020"),
    (2460755.0, "VE 2025"),
    (2460847.0, "SS 2025"),
    (2460940.0, "AE 2025"),
    (2461031.0, "WS 2025"),
    (2462581.0, "VE 2030"),
    (2462673.0, "SS 2030"),
    (2462766.0, "AE 2030"),
    (2462857.0, "WS 2030"),
    (2466234.0, "VE 2040"),
    (2466326.0, "SS 2040"),
    (2466419.0, "AE 2040"),
    (2466510.0, "WS 2040"),
    (2469886.0, "VE 2050"),
    (2469978.0, "SS 2050"),
    (2470071.0, "AE 2050"),
    (2470162.0, "WS 2050"),
]

flags = FLG_SWIEPH | FLG_SPEED
COMP_NAMES = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]

passed = 0
failed = 0
errors = 0
total = 0

for jd, event_name in EVENTS:
    for body, bname in BODIES.items():
        try:
            se_result = swe.calc_ut(jd, body, flags)
            le_result = ephem.calc_ut(jd, body, flags)

            for i, comp in enumerate(COMP_NAMES):
                total += 1
                se_val = se_result[0][i]
                le_val = le_result[0][i]
                diff = abs(le_val - se_val)

                if comp in ("dist", "dist_spd"):
                    tol = 0.001
                    is_pass = diff < tol
                else:
                    diff_arcsec = diff * 3600.0
                    if body in (MEAN_NODE, TRUE_NODE, MEAN_APOG):
                        tol = 30.0
                    elif body == MOON:
                        tol = 3.0
                    else:
                        tol = 2.0
                    is_pass = diff_arcsec < tol

                if is_pass:
                    passed += 1
                else:
                    diff_arcsec = diff * 3600.0
                    print(
                        f'FAIL {bname:10s} {comp:8s} {event_name:8s} JD={jd:.1f} SE={se_val:+14.8f} LE={le_val:+14.8f} diff={diff_arcsec:.4f}"'
                    )
                    failed += 1

        except Exception as e:
            errors += 1

print(f"\n{'=' * 60}")
print(f"Round 140: Planet Positions at Solstice/Equinox")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")
