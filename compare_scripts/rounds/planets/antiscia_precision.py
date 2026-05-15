"""Round 145: Planetary Antiscia Precision.

Tests antiscia (mirror points across 0°Cancer/0°Capricorn axis)
and contra-antiscia (mirror across 0°Aries/0°Libra axis) calculations.
Since these are derived from longitude, this validates longitude precision
from a different angle (pun intended).
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
)

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
ephem.set_ephe_path(_REF_EPHE_PATH)

BODIES = [
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
]
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
}

TEST_JDS = [
    2433282.5,
    2437000.5,
    2440587.5,
    2444239.5,
    2447892.5,
    2451545.0,
    2453371.5,
    2455197.5,
    2457023.5,
    2458849.5,
    2460310.5,
    2460676.5,
    2461041.5,
    2462502.5,
    2464328.5,
    2466154.5,
    2467980.5,
    2469807.5,
]

flags = FLG_SWIEPH | FLG_SPEED

passed = 0
failed = 0
total = 0

for jd in TEST_JDS:
    for body in BODIES:
        total += 1
        bname = BODY_NAMES[body]

        try:
            se_result = swe.calc_ut(jd, body, flags)
            le_result = ephem.calc_ut(jd, body, flags)

            se_lon = se_result[0][0]
            le_lon = le_result[0][0]

            # Antiscia: mirror across Cancer/Capricorn (0°-180° axis)
            # antiscia_lon = (360 - lon) % 360  -- simplified
            # Actually: antiscia = (180 - lon) % 360 for Cancer axis
            # No wait: antiscia across solstice axis = 360 - lon (reflected in 0°Cap/0°Can)
            # More precisely: the antiscion of X° is at (360° - X°) % 360°

            se_antiscia = (360.0 - se_lon) % 360.0
            le_antiscia = (360.0 - le_lon) % 360.0

            # Contra-antiscia: mirror across equinox axis (0°Ari/0°Lib)
            se_contra = (180.0 - se_lon) % 360.0
            le_contra = (180.0 - le_lon) % 360.0

            # The antiscia/contra-antiscia precision depends entirely on
            # the underlying longitude precision
            antiscia_diff = abs(le_antiscia - se_antiscia) * 3600.0
            contra_diff = abs(le_contra - se_contra) * 3600.0

            # Handle 360° wrap
            if antiscia_diff > 648000:  # > 180°
                antiscia_diff = 1296000 - antiscia_diff
            if contra_diff > 648000:
                contra_diff = 1296000 - contra_diff

            tol = 2.0  # 2" — same as longitude tolerance
            if body == MOON:
                tol = 3.0

            if antiscia_diff < tol and contra_diff < tol:
                passed += 1
            else:
                failed += 1
                print(
                    f'FAIL {bname:10s} JD={jd:.1f} antiscia_diff={antiscia_diff:.4f}" contra_diff={contra_diff:.4f}"'
                )

        except Exception as e:
            failed += 1
            print(f"ERR  {bname:10s} JD={jd:.1f}: {e}")

# Test 2: Verify antiscia pairs (planet A's antiscion should aspect planet B's position)
# This is a consistency check - if SE and LE agree on longitudes, their
# antiscia computations will also agree
print("\n=== Test 2: Antiscia pair consistency ===")

jd = 2460676.5  # 2025-01-01
all_se_lons = {}
all_le_lons = {}

for body in BODIES:
    total += 1
    try:
        se_result = swe.calc_ut(jd, body, flags)
        le_result = ephem.calc_ut(jd, body, flags)
        all_se_lons[body] = se_result[0][0]
        all_le_lons[body] = le_result[0][0]

        # Check each planet's antiscion
        se_anti = (360.0 - se_result[0][0]) % 360.0
        le_anti = (360.0 - le_result[0][0]) % 360.0

        diff = abs(le_anti - se_anti) * 3600.0
        if diff > 648000:
            diff = 1296000 - diff

        if diff < 3.0:
            passed += 1
        else:
            failed += 1
            print(f'FAIL pair {BODY_NAMES[body]:10s} diff={diff:.4f}"')

    except Exception as e:
        failed += 1

print(f"\n{'=' * 60}")
print(f"Round 145: Planetary Antiscia Precision")
print(f"{'=' * 60}")
print(f"Total:   {total}")
print(f"Passed:  {passed} ({100 * passed / max(total, 1):.1f}%)")
print(f"Failed:  {failed}")
