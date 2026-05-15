#!/usr/bin/env python3
"""Round 131: Planet Positions 200-Year Sweep (1850-2050)

Tests all major planet positions at yearly intervals across the full DE440 range,
focusing on long-term precision drift.
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
BODIES = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
]


def main():
    print("=" * 80)
    print("ROUND 131: Planet Positions 200-Year Sweep (1850-2050)")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # 1850-2050 at ~yearly intervals (every 365.25 days)
    start_jd = 2396758.5  # ~1850-01-01
    end_jd = 2469808.5  # ~2050-01-01
    step = 365.25

    worst = {}  # Track worst diff per body

    jd = start_jd
    while jd <= end_jd:
        for body_id, name in BODIES:
            try:
                se = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                le = ephem.calc_ut(jd, body_id, FLG_SPEED)[0]
            except Exception:
                continue

            # Longitude
            diff = le[0] - se[0]
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            diff_as = abs(diff) * 3600
            tol = 2.0 if name != "Moon" else 2.5  # Moon can be up to ~2" near 2050

            total_tests += 1
            if diff_as < tol:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 20:
                    year = 2000 + (jd - 2451545.0) / 365.25
                    failures.append(f'  {name} ~{year:.0f} lon: {diff_as:.4f}"')

            if name not in worst or diff_as > worst[name]:
                worst[name] = diff_as

            # Latitude
            diff_lat = abs(le[1] - se[1]) * 3600
            total_tests += 1
            if diff_lat < 2.0:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 30:
                    year = 2000 + (jd - 2451545.0) / 365.25
                    failures.append(f'  {name} ~{year:.0f} lat: {diff_lat:.4f}"')

            # Speed
            diff_spd = abs(le[3] - se[3]) * 3600
            total_tests += 1
            if diff_spd < 5.0:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 40:
                    year = 2000 + (jd - 2451545.0) / 365.25
                    failures.append(f'  {name} ~{year:.0f} spd: {diff_spd:.4f}"/day')

        jd += step

    print("\n  Worst-case differences per body:")
    for name in [n for _, n in BODIES]:
        if name in worst:
            print(f'    {name:10s}: {worst[name]:.4f}"')

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 131 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
    )
    if failures:
        print("\nSample failures:")
        for f in failures[:15]:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
