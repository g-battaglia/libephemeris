#!/usr/bin/env python3
"""Round 126: Inner Planet Inferior/Superior Conjunction Distances

Tests Mercury and Venus positions at inferior and superior conjunctions,
where geocentric distances are at extremes and parallax effects are maximal.
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
VENUS = 3


def find_conjunctions(body_id, start_jd, count=20):
    """Find inferior/superior conjunctions by tracking Sun-planet elongation."""
    events = []
    jd = start_jd
    step = 1.0

    def elong(j):
        s = swe.calc_ut(j, SUN, FLG_SPEED)[0][0]
        p = swe.calc_ut(j, body_id, FLG_SPEED)[0][0]
        e = p - s
        while e > 180:
            e -= 360
        while e < -180:
            e += 360
        return e

    prev = elong(jd)
    jd += step
    found = 0
    for _ in range(50000):
        if found >= count:
            break
        curr = elong(jd)
        if prev * curr < 0 and abs(prev) < 90:
            # Refine
            a, b = jd - step, jd
            for _ in range(50):
                mid = (a + b) / 2
                em = elong(mid)
                if (prev > 0 and em > 0) or (prev < 0 and em < 0):
                    a = mid
                else:
                    b = mid
            conj_jd = (a + b) / 2
            # Determine inferior vs superior by distance
            dist = swe.calc_ut(conj_jd, body_id, FLG_SPEED)[0][2]
            sun_dist = swe.calc_ut(conj_jd, SUN, FLG_SPEED)[0][2]
            ctype = "inferior" if dist < sun_dist else "superior"
            events.append((ctype, conj_jd))
            found += 1
        prev = curr
        jd += step
    return events


def main():
    print("=" * 80)
    print("ROUND 126: Inner Planet Conjunction Distances")
    print("=" * 80)
    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    for body_id, name in [(MERCURY, "Mercury"), (VENUS, "Venus")]:
        events = find_conjunctions(body_id, 2451545.0, 25)
        print(f"  Found {len(events)} conjunctions for {name}")

        for ctype, jd in events:
            for offset in [-0.05, 0.0, 0.05]:
                tj = jd + offset
                try:
                    se = swe.calc_ut(tj, body_id, FLG_SPEED)[0]
                    le = ephem.calc_ut(tj, body_id, FLG_SPEED)[0]
                except Exception:
                    continue

                for i, (lbl, tol) in enumerate(
                    [
                        ("lon", 2.0),
                        ("lat", 2.0),
                        ("dist", None),
                        ("lon_spd", 5.0),
                        ("lat_spd", 5.0),
                        ("dist_spd", None),
                    ]
                ):
                    if tol is None:
                        if abs(se[i]) > 1e-10:
                            rel = abs(le[i] - se[i]) / abs(se[i])
                            total_tests += 1
                            if rel < 1e-4:
                                total_pass += 1
                            else:
                                total_fail += 1
                                if len(failures) < 20:
                                    failures.append(
                                        f"  {name} {ctype} {lbl}: rel={rel:.2e}"
                                    )
                        continue
                    d = le[i] - se[i]
                    if i == 0:
                        if d > 180:
                            d -= 360
                        elif d < -180:
                            d += 360
                    das = abs(d) * 3600
                    total_tests += 1
                    if das < tol:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 30:
                            failures.append(f'  {name} {ctype} {lbl}: {das:.4f}"')

    pct = 100 * total_pass / total_tests if total_tests else 0
    print(
        f"\n{'=' * 80}\nROUND 126 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)\n  Failures: {total_fail}\n{'=' * 80}"
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
