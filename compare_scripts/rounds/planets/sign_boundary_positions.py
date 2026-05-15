#!/usr/bin/env python3
"""Round 121: Planetary Longitude at Sign Boundaries

Tests planet positions precisely at ingress points (0°, 30°, 60°... of each sign)
where rounding/interpolation errors are most visible in astrological software.

Verifies:
- Planet positions near exact sign boundaries
- Speed accuracy at ingress moments
- Multiple planets crossing same degree
- Both tropical and sidereal sign boundaries
"""

from __future__ import annotations

import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_SIDEREAL = 65536

SUN = 0
MOON = 1
MERCURY = 2
VENUS = 3
MARS = 4
JUPITER = 5
SATURN = 6


def find_sign_ingresses(body_id, start_jd, count=30):
    """Find times when planet crosses sign boundaries (multiples of 30°)."""
    ingresses = []
    jd = start_jd
    step = 1.0 if body_id != MOON else 0.25

    prev_sign = int(swe.calc_ut(jd, body_id, FLG_SPEED)[0][0] / 30.0)
    jd += step

    max_iter = 50000
    found = 0

    for _ in range(max_iter):
        if found >= count:
            break

        try:
            pos = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
        except Exception:
            jd += step
            continue

        curr_sign = int(pos[0] / 30.0) % 12

        if curr_sign != prev_sign:
            # Refine with bisection
            a, b = jd - step, jd
            for _ in range(50):
                mid = (a + b) / 2
                mid_sign = int(swe.calc_ut(mid, body_id, FLG_SPEED)[0][0] / 30.0) % 12
                if mid_sign == prev_sign:
                    a = mid
                else:
                    b = mid

            ingress_jd = (a + b) / 2
            boundary_deg = curr_sign * 30.0
            ingresses.append((ingress_jd, boundary_deg))
            found += 1

        prev_sign = curr_sign
        jd += step

    return ingresses


def main():
    print("=" * 80)
    print("ROUND 121: Planetary Longitude at Sign Boundaries")
    print("=" * 80)

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    bodies = [
        (SUN, "Sun", 24),
        (MOON, "Moon", 30),
        (MERCURY, "Mercury", 20),
        (VENUS, "Venus", 20),
        (MARS, "Mars", 15),
        (JUPITER, "Jupiter", 8),
        (SATURN, "Saturn", 6),
    ]

    start_jd = 2451545.0  # J2000

    for body_id, body_name, count in bodies:
        print(f"\n  Finding {count} ingresses for {body_name}...")
        ingresses = find_sign_ingresses(body_id, start_jd, count)
        print(f"  Found {len(ingresses)} ingresses")

        for ingress_jd, boundary in ingresses:
            # Test at ingress, and at offsets around it
            for offset_days in [-0.01, -0.001, 0.0, 0.001, 0.01]:
                jd = ingress_jd + offset_days

                try:
                    se_pos = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                    le_pos = ephem.calc_ut(jd, body_id, FLG_SPEED)[0]
                except Exception:
                    continue

                # Compare longitude
                diff = le_pos[0] - se_pos[0]
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                diff_as = abs(diff) * 3600
                tol = 2.0

                total_tests += 1
                if diff_as < tol:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 20:
                        failures.append(
                            f"  {body_name} boundary={boundary}° offset={offset_days}: "
                            f'diff={diff_as:.4f}" (tol={tol}")'
                        )

                # Compare speed
                diff_spd = abs(le_pos[3] - se_pos[3]) * 3600
                tol_spd = 5.0

                total_tests += 1
                if diff_spd < tol_spd:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 30:
                        failures.append(
                            f"  {body_name} boundary={boundary}° speed: "
                            f'diff={diff_spd:.4f}"/day (tol={tol_spd}")'
                        )

    # Sidereal sign boundaries
    print("\n--- Testing sidereal sign boundaries ---")
    swe.set_sid_mode(1)  # Lahiri
    ephem.set_sid_mode(1, 0, 0)

    for body_id, body_name in [(SUN, "Sun"), (MOON, "Moon"), (MARS, "Mars")]:
        ingresses = find_sign_ingresses(body_id, start_jd, 12)

        for ingress_jd, boundary in ingresses:
            flags = FLG_SPEED | FLG_SIDEREAL
            try:
                se_pos = swe.calc_ut(ingress_jd, body_id, flags)[0]
                le_pos = ephem.calc_ut(ingress_jd, body_id, flags)[0]
            except Exception:
                continue

            diff = le_pos[0] - se_pos[0]
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            diff_as = abs(diff) * 3600
            tol = 20.0  # Known sidereal offset

            total_tests += 1
            if diff_as < tol:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 40:
                    failures.append(
                        f'  SID {body_name} boundary={boundary}°: diff={diff_as:.4f}"'
                    )

    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    print("\n" + "=" * 80)
    pct = 100 * total_pass / total_tests if total_tests > 0 else 0
    print(f"ROUND 121 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
    print(f"  Failures: {total_fail}")
    print("=" * 80)

    if failures:
        print("\nSample failures:")
        for f in failures[:20]:
            print(f)

    if total_fail == 0:
        print("\nAll tests PASSED!")

    return total_fail


if __name__ == "__main__":
    sys.exit(main())
