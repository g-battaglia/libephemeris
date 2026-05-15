"""Round 135: Planetary Latitude Speed Accuracy Deep Sweep.

Tests latitude speed (dlat/dt) across all planets, multiple epochs,
and flag combinations. Latitude speed is often the weakest link due
to finite-difference methodology differences.
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
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
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

BODY_NAMES = {
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
    MEAN_NODE: "MeanNode",
    TRUE_NODE: "TrueNode",
    MEAN_APOG: "MeanApog",
    OSCU_APOG: "OscuApog",
    CHIRON: "Chiron",
    CERES: "Ceres",
    PALLAS: "Pallas",
    JUNO: "Juno",
    VESTA: "Vesta",
}

# Test dates spanning 200 years
TEST_JDS = [
    2378496.5,  # 1800-01-01
    2396758.5,  # 1850-01-01
    2415020.5,  # 1900-01-01
    2433282.5,  # 1950-01-01
    2440587.5,  # 1970-01-01
    2444239.5,  # 1980-01-01
    2447892.5,  # 1990-01-01
    2451545.0,  # 2000-01-01.5 (J2000)
    2453371.5,  # 2005-01-01
    2455197.5,  # 2010-01-01
    2457023.5,  # 2015-01-01
    2458849.5,  # 2020-01-01
    2460310.5,  # 2024-01-01
    2460676.5,  # 2025-01-01
    2462502.5,  # 2030-01-01
    2469807.5,  # 2050-01-01
    2477113.5,  # 2070-01-01
    2488069.5,  # 2100-01-01
]

# Flag combinations
FLAG_COMBOS = {
    "default": FLG_SWIEPH | FLG_SPEED,
    "J2000": FLG_SWIEPH | FLG_SPEED | FLG_J2000,
    "NONUT": FLG_SWIEPH | FLG_SPEED | FLG_NONUT,
    "TRUEPOS": FLG_SWIEPH | FLG_SPEED | FLG_TRUEPOS,
    "NOABERR": FLG_SWIEPH | FLG_SPEED | FLG_NOABERR,
    "HELCTR": FLG_SWIEPH | FLG_SPEED | FLG_HELCTR,
}

# Bodies that can't be heliocentric
NO_HELIO = {SUN, MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG}

# Bodies not available at all dates
RESTRICTED_BODIES = {CHIRON, CERES, PALLAS, JUNO, VESTA}
RESTRICTED_MIN_JD = 2396758.5  # ~1850

passed = 0
failed = 0
errors = 0
total = 0

worst_cases = []

for jd in TEST_JDS:
    for body, bname in BODY_NAMES.items():
        if body in RESTRICTED_BODIES and jd < RESTRICTED_MIN_JD:
            continue

        for flag_name, flags in FLAG_COMBOS.items():
            if flag_name == "HELCTR" and body in NO_HELIO:
                continue
            # Sun can't be heliocentric
            if (flags & FLG_HELCTR) and body == SUN:
                continue

            total += 1
            try:
                se_result = swe.calc_ut(jd, body, flags)
                se_lat_speed = se_result[0][4]  # lat_speed is index 4

                le_result = ephem.calc_ut(jd, body, flags)
                le_lat_speed = le_result[0][4]

                diff = abs(le_lat_speed - se_lat_speed)
                # Convert to arcseconds
                diff_arcsec = diff * 3600.0

                # Adaptive tolerance based on body type and speed magnitude
                base_tol = 1.0  # 1 arcsec/day default

                # Analytical bodies have known differences
                if body in (MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG):
                    base_tol = 60.0  # 1 arcmin/day
                elif body == MOON:
                    base_tol = 5.0  # Moon moves fast, finite-diff differences
                elif body in (CHIRON,):
                    base_tol = 2.0

                # Relative tolerance: 0.5% of speed magnitude
                speed_mag = max(abs(se_lat_speed), abs(le_lat_speed))
                rel_tol = speed_mag * 3600.0 * 0.005  # 0.5%

                tol = max(base_tol, rel_tol)

                if diff_arcsec <= tol:
                    passed += 1
                else:
                    failed += 1
                    info = f'FAIL {bname:10s} JD={jd} {flag_name:10s} SE={se_lat_speed:+12.6f} LE={le_lat_speed:+12.6f} diff={diff_arcsec:.3f}"/day tol={tol:.3f}"/day'
                    print(info)
                    worst_cases.append((diff_arcsec, info))

            except Exception as e:
                errors += 1
                print(f"ERR  {bname:10s} JD={jd} {flag_name:10s}: {e}")

# Summary
print(f"\n{'=' * 60}")
print(f"Round 135: Planetary Latitude Speed Accuracy")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
print(f"Errors:  {errors}")

if worst_cases:
    worst_cases.sort(key=lambda x: -x[0])
    print(f"\nTop 10 worst cases:")
    for diff, info in worst_cases[:10]:
        print(f"  {info}")
