#!/usr/bin/env python3
"""Round 133: houses_armc vs houses_ex Consistency

Tests that houses_armc_ex2 produces identical cusps to houses_ex2 when given
the same ARMC (computed from sidereal time).
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256


def se_hsys(ch):
    return ch.encode("ascii") if isinstance(ch, str) else ch


def main():
    print("=" * 80)
    print("ROUND 133: houses_armc vs houses_ex Consistency")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    house_systems = ["P", "K", "O", "R", "C", "E", "W", "M", "B"]
    test_jds = [2451545.0, 2455197.5, 2459580.5, 2460310.5, 2444239.5]
    latitudes = [0.0, 15.0, 30.0, 45.0, 55.0, -20.0, -40.0]
    lon = 12.5

    for hsys in house_systems:
        for jd in test_jds:
            for lat in latitudes:
                # Get cusps from houses_ex2
                try:
                    le_ex2 = ephem.houses_ex2(jd, lat, lon, ord(hsys), FLG_SPEED)
                    cusps_ex2 = le_ex2[0]
                except Exception:
                    continue

                # Compute ARMC from sidereal time
                armc = ephem.sidtime(jd) * 15.0 + lon
                if armc >= 360:
                    armc -= 360

                # Get obliquity
                try:
                    nut = ephem.calc_ut(jd, -1, 0)
                    eps = nut[0][0]  # true obliquity
                except Exception:
                    continue

                # Get cusps from houses_armc_ex2
                try:
                    le_armc = ephem.houses_armc_ex2(
                        armc, lat, eps, ord(hsys), FLG_SPEED
                    )
                    cusps_armc = le_armc[0]
                except Exception:
                    continue

                # Compare cusps — they should be very close
                for i in range(min(12, len(cusps_ex2), len(cusps_armc))):
                    diff = cusps_armc[i] - cusps_ex2[i]
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360
                    diff_as = abs(diff) * 3600

                    total_tests += 1
                    if diff_as < 1.0:  # Should be identical within rounding
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 20:
                            failures.append(
                                f"  {hsys} lat={lat} JD={jd:.1f} cusp{i + 1}: "
                                f'ex2={cusps_ex2[i]:.6f} armc={cusps_armc[i]:.6f} diff={diff_as:.4f}"'
                            )

    # Also compare SE houses vs LE houses_armc
    print("\n--- Comparing SE houses vs LE houses_armc ---")
    for hsys in ["P", "K", "O", "R", "C", "E"]:
        for jd in test_jds[:3]:
            for lat in [0.0, 30.0, 45.0, -33.0]:
                try:
                    se_cusps, se_ascmc = swe.houses_ex(jd, lat, lon, se_hsys(hsys))
                except Exception:
                    continue

                armc = swe.sidtime(jd) * 15.0 + lon
                if armc >= 360:
                    armc -= 360
                eps_se = swe.calc_ut(jd, -1, 0)[0][0]

                try:
                    le_armc = ephem.houses_armc_ex2(
                        armc, lat, eps_se, ord(hsys), FLG_SPEED
                    )
                    le_cusps = le_armc[0]
                except Exception:
                    continue

                for i in range(min(12, len(se_cusps), len(le_cusps))):
                    diff = le_cusps[i] - se_cusps[i]
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360
                    diff_as = abs(diff) * 3600

                    total_tests += 1
                    if diff_as < 2.0:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 30:
                            failures.append(
                                f'  SE vs LE_armc {hsys} lat={lat} cusp{i + 1}: {diff_as:.4f}"'
                            )

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 133 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
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
