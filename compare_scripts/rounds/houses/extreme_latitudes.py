#!/usr/bin/env python3
"""Round 96: House Systems at Extreme Latitudes — Deep Sweep

Tests house cusp calculations at challenging latitudes:
1. Arctic/Antarctic circles (66.5°)
2. Near-polar (70°, 75°, 80°, 85°)
3. Equator (0°)
4. Tropical latitudes (23.4°)
5. Multiple house systems (P, K, O, R, C, E, W, X, M, B, G, H, T, U, Y)
6. Multiple ARMC values (0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)
7. Both hemispheres
"""

from __future__ import annotations

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

# House systems to test
HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal"),
    ("W", "Whole Sign"),
    ("X", "Meridian"),
    ("M", "Morinus"),
    ("B", "Alcabitius"),
    ("T", "Polich-Page"),
    ("U", "Krusinski"),
    ("H", "Horizon"),
    ("G", "Gauquelin"),
]

# Latitudes to test (both hemispheres)
LATITUDES = [
    0.0,
    23.4,
    45.0,
    55.0,
    60.0,
    66.0,
    66.5,
    67.0,
    70.0,
    75.0,
    80.0,
    85.0,
    89.0,
]

# ARMC values
ARMC_VALUES = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]

# Obliquity at J2000
OBLIQUITY = 23.4392911

# Tolerances
CUSP_TOL = 0.01  # degrees (36") for cusps
ASCMC_TOL = 0.01  # degrees for ASC/MC/etc
KNOWN_ISSUES = {
    "I": True,  # Sunshine house system - known differences
}


def se_hsys(ch):
    return ch.encode("ascii")


def le_hsys(ch):
    return ord(ch)


def run_tests():
    passed = 0
    failed = 0
    errors = 0
    total = 0
    fail_details = []

    print("=" * 80)
    print("ROUND 96: House Systems at Extreme Latitudes — Deep Sweep")
    print("=" * 80)

    for hsys_ch, hsys_name in HOUSE_SYSTEMS:
        h_pass = 0
        h_fail = 0
        h_err = 0

        for lat_abs in LATITUDES:
            for lat_sign in [1, -1]:
                lat = lat_abs * lat_sign
                if lat == 0.0 and lat_sign == -1:
                    continue  # Skip duplicate equator

                for armc in ARMC_VALUES:
                    total += 1

                    try:
                        se_cusps, se_ascmc = swe.houses_armc(
                            armc, lat, OBLIQUITY, se_hsys(hsys_ch)
                        )
                    except Exception as e:
                        h_err += 1
                        errors += 1
                        continue

                    try:
                        le_result = ephem.houses_armc(
                            armc, lat, OBLIQUITY, le_hsys(hsys_ch)
                        )
                        le_cusps = le_result[0]
                        le_ascmc = le_result[1]
                    except Exception as e:
                        h_err += 1
                        errors += 1
                        if h_err <= 3:
                            fail_details.append(
                                f"  ERROR [{hsys_name}] lat={lat} armc={armc}: {e}"
                            )
                        continue

                    # Compare cusps
                    n_cusps = min(len(se_cusps), len(le_cusps))
                    cusp_ok = True
                    for i in range(n_cusps):
                        if se_cusps[i] == 0.0 and le_cusps[i] == 0.0:
                            continue
                        diff = abs(se_cusps[i] - le_cusps[i])
                        if diff > 180:
                            diff = 360 - diff
                        if diff > CUSP_TOL:
                            cusp_ok = False
                            if len(fail_details) < 20:
                                fail_details.append(
                                    f"  FAIL [{hsys_name}] lat={lat:6.1f} armc={armc:3d} "
                                    f"cusp[{i + 1}]: SE={se_cusps[i]:.6f} LE={le_cusps[i]:.6f} "
                                    f'diff={diff * 3600:.1f}"'
                                )
                            break

                    # Compare ASC, MC, ARMC, Vertex (indices 0-3)
                    ascmc_ok = True
                    ascmc_names = ["ASC", "MC", "ARMC", "Vertex"]
                    for i in range(min(4, len(se_ascmc), len(le_ascmc))):
                        if se_ascmc[i] == 0.0 and le_ascmc[i] == 0.0:
                            continue
                        diff = abs(se_ascmc[i] - le_ascmc[i])
                        if diff > 180:
                            diff = 360 - diff
                        if diff > ASCMC_TOL:
                            ascmc_ok = False
                            if len(fail_details) < 20:
                                fail_details.append(
                                    f"  FAIL [{hsys_name}] lat={lat:6.1f} armc={armc:3d} "
                                    f"{ascmc_names[i]}: SE={se_ascmc[i]:.6f} "
                                    f'LE={le_ascmc[i]:.6f} diff={diff * 3600:.1f}"'
                                )
                            break

                    if cusp_ok and ascmc_ok:
                        h_pass += 1
                        passed += 1
                    else:
                        h_fail += 1
                        failed += 1

        pct = 100 * h_pass / max(1, h_pass + h_fail + h_err)
        status = "OK" if h_fail == 0 else "ISSUES"
        print(
            f"  [{hsys_name:15s}] {h_pass}/{h_pass + h_fail + h_err} passed ({pct:.1f}%) "
            f"[{h_fail} fail, {h_err} err] {status}"
        )

    print()
    if fail_details:
        print("Sample failures:")
        for d in fail_details[:20]:
            print(d)

    print()
    print("=" * 80)
    print(
        f"ROUND 96 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%), "
        f"{failed} failed, {errors} errors"
    )
    print("=" * 80)


if __name__ == "__main__":
    t0 = time.time()
    run_tests()
    print(f"\nCompleted in {time.time() - t0:.1f}s")
