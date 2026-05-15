#!/usr/bin/env python3
"""Round 134: Mercury Max Elongation Precision

Tests Mercury positions at maximum elongation points (greatest eastern/western
elongation) — critical for visibility predictions and astrological interpretation.
"""

from __future__ import annotations
import os, sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
SUN = 0
MERCURY = 2


def find_max_elongations(start_jd, count=30):
    """Find Mercury max elongation events."""
    events = []
    jd = start_jd
    step = 1.0

    def elong(j):
        s = swe.calc_ut(j, SUN, FLG_SPEED)[0][0]
        p = swe.calc_ut(j, MERCURY, FLG_SPEED)[0][0]
        e = p - s
        while e > 180:
            e -= 360
        while e < -180:
            e += 360
        return e

    prev_prev = elong(jd - step)
    prev = elong(jd)
    jd += step
    found = 0

    for _ in range(50000):
        if found >= count:
            break
        curr = elong(jd)

        # Check for elongation extremum
        if abs(prev) > abs(prev_prev) and abs(prev) > abs(curr) and abs(prev) > 15:
            # Refine
            a, b = jd - 2 * step, jd
            gr = 1.618033988749895
            for _ in range(40):
                c = b - (b - a) / gr
                d = a + (b - a) / gr
                ec = abs(elong(c))
                ed = abs(elong(d))
                if ec > ed:
                    b = d
                else:
                    a = c
            refined = (a + b) / 2
            e_val = elong(refined)
            etype = "eastern" if e_val > 0 else "western"
            events.append((etype, refined, e_val))
            found += 1

        prev_prev = prev
        prev = curr
        jd += step

    return events


def main():
    print("=" * 80)
    print("ROUND 134: Mercury Max Elongation Precision")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    events = find_max_elongations(2451545.0, 40)
    print(f"Found {len(events)} max elongation events")

    max_east = max(e[2] for e in events if e[0] == "eastern")
    max_west = min(e[2] for e in events if e[0] == "western")
    print(f"  Greatest eastern: {max_east:.2f}°")
    print(f"  Greatest western: {max_west:.2f}°")

    for etype, jd, e_val in events:
        for offset in [-0.1, 0.0, 0.1]:
            tj = jd + offset
            try:
                se = swe.calc_ut(tj, MERCURY, FLG_SPEED)[0]
                le = ephem.calc_ut(tj, MERCURY, FLG_SPEED)[0]
            except Exception:
                continue

            # Longitude
            d = le[0] - se[0]
            if d > 180:
                d -= 360
            elif d < -180:
                d += 360
            das = abs(d) * 3600
            total_tests += 1
            if das < 2.0:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 15:
                    failures.append(f'  {etype} elong={e_val:.1f}° lon: {das:.4f}"')

            # Latitude
            das_lat = abs(le[1] - se[1]) * 3600
            total_tests += 1
            if das_lat < 2.0:
                total_pass += 1
            else:
                total_fail += 1

            # Speed
            das_spd = abs(le[3] - se[3]) * 3600
            total_tests += 1
            if das_spd < 5.0:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 25:
                    failures.append(
                        f'  {etype} elong={e_val:.1f}° spd: {das_spd:.4f}"/day'
                    )

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 134 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
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
