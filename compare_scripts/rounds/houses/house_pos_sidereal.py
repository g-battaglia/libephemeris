"""Round 142: house_pos with Sidereal Mode.

Tests house_pos() with sidereal mode enabled across multiple
house systems, locations, and ayanamsha modes.
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
    FLG_SIDEREAL,
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
)

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

BODIES = [SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN]
BODY_NAMES = {
    SUN: "Sun",
    MOON: "Moon",
    MERCURY: "Mercury",
    VENUS: "Venus",
    MARS: "Mars",
    JUPITER: "Jupiter",
    SATURN: "Saturn",
}

TEST_JDS = [2451545.0, 2455197.5, 2458849.5, 2460676.5, 2462502.5]

LOCATIONS = [
    (51.5, -0.1, "London"),
    (40.7, -74.0, "NewYork"),
    (28.6, 77.2, "Delhi"),
    (-33.9, 151.2, "Sydney"),
    (35.7, 139.7, "Tokyo"),
    (0.0, 0.0, "Equator"),
    (64.1, -21.9, "Reykjavik"),
]

HOUSE_SYSTEMS = [b"P", b"K", b"O", b"R", b"C", b"E", b"W", b"B"]
LE_HOUSE_SYSTEMS = [ord(c) for c in ["P", "K", "O", "R", "C", "E", "W", "B"]]

SIDEREAL_MODES = [0, 1, 3, 27]  # Fagan, Lahiri, Raman, TrueCitra

passed = 0
failed = 0
errors = 0
total = 0

for sid_mode in SIDEREAL_MODES:
    swe.set_sid_mode(sid_mode)
    ephem.set_sid_mode(sid_mode)

    for jd in TEST_JDS:
        # Get obliquity
        flags = FLG_SWIEPH | FLG_SPEED | FLG_SIDEREAL

        for lat, lon, loc_name in LOCATIONS:
            # Compute ARMC for this location and time
            try:
                se_sidtime = swe.sidtime(jd)
                armc = (se_sidtime * 15.0 + lon) % 360.0
                eps = swe.calc_ut(jd, -1, FLG_SWIEPH)[0][
                    0
                ]  # obliquity via special body
            except Exception:
                # Fallback: compute obliquity
                eps = 23.4393

            for body in BODIES:
                # Get planet position in sidereal
                try:
                    se_pos = swe.calc_ut(jd, body, flags)
                    le_pos = ephem.calc_ut(jd, body, flags)
                    plon = se_pos[0][0]
                    plat = se_pos[0][1]
                    le_plon = le_pos[0][0]
                    le_plat = le_pos[0][1]
                except Exception:
                    continue

                for hi, (se_hsys, le_hsys) in enumerate(
                    zip(HOUSE_SYSTEMS, LE_HOUSE_SYSTEMS)
                ):
                    total += 1
                    try:
                        se_hp = swe.house_pos(armc, lat, eps, (plon, plat), se_hsys)
                        le_hp = ephem.house_pos(
                            armc, lat, eps, le_hsys, le_plon, le_plat
                        )

                        diff = abs(se_hp - le_hp)
                        # Wrap around 12 houses
                        if diff > 6.0:
                            diff = 12.0 - diff

                        tol = 0.01  # 0.01 house units (~0.3°)

                        if diff < tol:
                            passed += 1
                        else:
                            failed += 1
                            hsys_ch = chr(le_hsys)
                            print(
                                f"FAIL {BODY_NAMES[body]:8s} {hsys_ch} sid={sid_mode} {loc_name:10s} JD={jd:.1f} SE={se_hp:.6f} LE={le_hp:.6f} diff={diff:.6f}"
                            )

                    except Exception as e:
                        errors += 1
                        if "polar" not in str(e).lower():
                            hsys_ch = chr(le_hsys)
                            print(
                                f"ERR  {BODY_NAMES[body]:8s} {hsys_ch} sid={sid_mode} {loc_name:10s}: {e}"
                            )

swe.set_sid_mode(0)
ephem.set_sid_mode(0)

print(f"\n{'=' * 60}")
print(f"Round 142: house_pos with Sidereal Mode")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")
