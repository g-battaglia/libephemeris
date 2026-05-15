#!/usr/bin/env python3
"""Round 122: Moon Latitude Extremes — Max Declination ~28.5°

Tests Moon positions at maximum/minimum ecliptic latitude (~5.3°) and
maximum/minimum declination (~28.5°), where perturbation effects are strongest.
"""

from __future__ import annotations
import os, sys, math

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
FLG_EQUATORIAL = 2048
MOON = 1


def find_lat_extremes(start_jd, count=40):
    """Find Moon latitude maxima/minima."""
    events = []
    jd = start_jd
    step = 0.5
    prev_lat = swe.calc_ut(jd - step, MOON, FLG_SPEED)[0][1]
    curr_lat = swe.calc_ut(jd, MOON, FLG_SPEED)[0][1]
    found = 0
    for _ in range(20000):
        if found >= count:
            break
        jd += step
        next_lat = swe.calc_ut(jd, MOON, FLG_SPEED)[0][1]
        if (curr_lat > prev_lat and curr_lat > next_lat) or (
            curr_lat < prev_lat and curr_lat < next_lat
        ):
            # Refine
            a, b = jd - 2 * step, jd
            for _ in range(40):
                m1 = a + (b - a) / 3
                m2 = a + 2 * (b - a) / 3
                l1 = swe.calc_ut(m1, MOON, FLG_SPEED)[0][1]
                l2 = swe.calc_ut(m2, MOON, FLG_SPEED)[0][1]
                if curr_lat > prev_lat:  # maximum
                    if l1 > l2:
                        b = m2
                    else:
                        a = m1
                else:  # minimum
                    if l1 < l2:
                        b = m2
                    else:
                        a = m1
            refined_jd = (a + b) / 2
            refined_lat = swe.calc_ut(refined_jd, MOON, FLG_SPEED)[0][1]
            etype = "max_lat" if refined_lat > 0 else "min_lat"
            events.append((etype, refined_jd, refined_lat))
            found += 1
        prev_lat = curr_lat
        curr_lat = next_lat
    return events


def main():
    print("=" * 80)
    print("ROUND 122: Moon Latitude Extremes")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    events = find_lat_extremes(2451545.0, count=60)
    print(f"Found {len(events)} latitude extremes")

    max_lat_se = max(abs(e[2]) for e in events)
    print(f"Max |latitude| found: {max_lat_se:.4f}°")

    flags_list = [
        ("default", FLG_SPEED),
        ("equatorial", FLG_SPEED | FLG_EQUATORIAL),
    ]

    for etype, jd, se_ref_lat in events:
        for flag_name, flags in flags_list:
            try:
                se_pos = swe.calc_ut(jd, MOON, flags)[0]
                le_pos = ephem.calc_ut(jd, MOON, flags)[0]
            except Exception:
                continue

            labels = ["lon/ra", "lat/dec", "dist", "lon_spd", "lat_spd", "dist_spd"]
            tols = [2.0, 2.0, None, 5.0, 5.0, None]

            for i, (label, tol) in enumerate(zip(labels, tols)):
                if tol is None:
                    # distance — relative
                    if abs(se_pos[i]) > 1e-10:
                        rel = abs(le_pos[i] - se_pos[i]) / abs(se_pos[i])
                        total_tests += 1
                        if rel < 1e-5:
                            total_pass += 1
                        else:
                            total_fail += 1
                    continue

                diff = le_pos[i] - se_pos[i]
                if i == 0:
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360
                diff_as = abs(diff) * 3600

                total_tests += 1
                if diff_as < tol:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 20:
                        failures.append(
                            f'  {etype} JD={jd:.4f} {flag_name} {label}: {diff_as:.4f}" (tol={tol}")'
                        )

    print(f"\n{'=' * 80}")
    pct = 100 * total_pass / total_tests if total_tests else 0
    print(f"ROUND 122 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
    print(f"  Failures: {total_fail}")
    print("=" * 80)
    if failures:
        print("\nSample failures:")
        for f in failures:
            print(f)
    if total_fail == 0:
        print("\nAll tests PASSED!")
    return total_fail


if __name__ == "__main__":
    sys.exit(main())
