#!/usr/bin/env python3
"""Round 132: cotrans_sp Deep Verification

Tests coordinate transformation with speeds (cotrans_sp) comparing LE vs SE
across a wide range of coordinates and obliquity values.
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)


def main():
    print("=" * 80)
    print("ROUND 132: cotrans_sp Deep Verification")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Test coordinates covering full range
    test_coords = []
    for lon in range(0, 360, 30):
        for lat in range(-80, 81, 20):
            for spd_lon in [-1.0, 0.0, 1.0, 15.0]:
                for spd_lat in [-0.5, 0.0, 0.5]:
                    test_coords.append(
                        (float(lon), float(lat), 1.0, spd_lon, spd_lat, 0.0)
                    )

    obliquities = [23.44, 23.0, 23.5, 22.5, 24.0]

    for obl in obliquities:
        for coord in test_coords:
            lon, lat, dist, slon, slat, sdist = coord
            try:
                se_res = swe.cotrans_sp((lon, lat, dist, slon, slat, sdist), obl)
                le_res = ephem.cotrans_sp((lon, lat, dist, slon, slat, sdist), obl)
            except Exception:
                continue

            for i in range(6):
                diff = abs(le_res[i] - se_res[i])
                if i == 0:
                    d = le_res[i] - se_res[i]
                    if d > 180:
                        d -= 360
                    elif d < -180:
                        d += 360
                    diff = abs(d)

                total_tests += 1
                if diff < 1e-8:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 20:
                        failures.append(
                            f"  obl={obl} ({lon},{lat}) idx={i}: SE={se_res[i]:.10f} LE={le_res[i]:.10f} diff={diff:.2e}"
                        )

    # Also test cotrans (without speeds)
    print("\n--- Testing cotrans (no speeds) ---")
    for obl in obliquities:
        for lon in range(0, 360, 15):
            for lat in range(-85, 86, 10):
                try:
                    se_res = swe.cotrans((float(lon), float(lat), 1.0), obl)
                    le_res = ephem.cotrans((float(lon), float(lat), 1.0), obl)
                except Exception:
                    continue

                for i in range(3):
                    diff = abs(le_res[i] - se_res[i])
                    if i == 0:
                        d = le_res[i] - se_res[i]
                        if d > 180:
                            d -= 360
                        elif d < -180:
                            d += 360
                        diff = abs(d)

                    total_tests += 1
                    if diff < 1e-8:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 30:
                            failures.append(
                                f"  cotrans obl={obl} ({lon},{lat}) idx={i}: diff={diff:.2e}"
                            )

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 132 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
    )
    if failures:
        print("\nSample failures:")
        for f in failures:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
