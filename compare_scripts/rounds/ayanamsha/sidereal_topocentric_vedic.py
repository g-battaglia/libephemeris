#!/usr/bin/env python3
"""Round 205: Sidereal+Topocentric combo.

Tests planet positions with both FLG_SIDEREAL and FLG_TOPOCTR flags
set simultaneously — a common real-world use case for Vedic astrology.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

BODIES = [
    ("Sun", ephem.SUN, swe.SUN),
    ("Moon", ephem.MOON, swe.MOON),
    ("Mars", ephem.MARS, swe.MARS),
    ("Jupiter", ephem.JUPITER, swe.JUPITER),
    ("Saturn", ephem.SATURN, swe.SATURN),
    ("MeanNode", ephem.MEAN_NODE, swe.MEAN_NODE),
]

LOCATIONS = [
    ("Delhi", 28.6139, 77.2090, 216.0),
    ("Mumbai", 19.0760, 72.8777, 14.0),
    ("London", 51.5074, -0.1278, 11.0),
    ("Sydney", -33.8688, 151.2093, 58.0),
]

SID_MODES = [
    ("Lahiri", 1),
    ("Raman", 3),
    ("Krishnamurti", 5),
]

TEST_JDS = [2451545.0, 2455197.5, 2458849.5, 2460310.5]


def test_sidereal_topo():
    global passed, failed, total

    print("=" * 70)
    print("Round 205: Sidereal + Topocentric Combo")
    print("=" * 70)

    for sid_name, sid_mode in SID_MODES:
        print(f"\n--- {sid_name} (mode {sid_mode}) ---")
        ephem.set_sid_mode(sid_mode, 0.0, 0.0)
        swe.set_sid_mode(sid_mode, 0.0, 0.0)

        for loc_name, lat, lon, alt in LOCATIONS:
            ephem.set_topo(lon, lat, alt)
            swe.set_topo(lon, lat, alt)

            flags_le = (
                ephem.FLG_SWIEPH
                | ephem.FLG_SPEED
                | ephem.FLG_SIDEREAL
                | ephem.FLG_TOPOCTR
            )
            flags_se = (
                swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_SIDEREAL | swe.FLG_TOPOCTR
            )

            for jd in TEST_JDS:
                for bname, le_b, se_b in BODIES:
                    try:
                        le_r = ephem.calc_ut(jd, le_b, flags_le)
                        se_r = swe.calc_ut(jd, se_b, flags_se)
                    except Exception:
                        continue

                    total += 1
                    lon_diff = abs(le_r[0][0] - se_r[0][0])
                    if lon_diff > 180:
                        lon_diff = 360 - lon_diff
                    lon_as = lon_diff * 3600

                    # Topocentric+Sidereal: ~15-20" tolerance (sidereal ~14" + topo parallax)
                    tol = 20.0 if bname == "Moon" else 15.0
                    if lon_as <= tol:
                        passed += 1
                    else:
                        failed += 1
                        failures.append(
                            f'  {sid_name} {loc_name} {bname} JD={jd:.1f}: diff={lon_as:.2f}"'
                        )

    # Reset
    ephem.set_sid_mode(0, 0.0, 0.0)
    swe.set_sid_mode(0, 0.0, 0.0)


if __name__ == "__main__":
    test_sidereal_topo()

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
