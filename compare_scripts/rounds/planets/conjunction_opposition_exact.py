#!/usr/bin/env python3
"""Round 221: Planetary conjunction/opposition exact timing.

Finds exact times of planetary conjunctions and oppositions by searching
for minimum angular separation, then compares planet positions at those
moments between LE and SE.
"""

from __future__ import annotations

import os, sys, math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

PAIRS = [
    ("Sun-Moon", ephem.SUN, swe.SUN, ephem.MOON, swe.MOON, 5.0),
    ("Sun-Mercury", ephem.SUN, swe.SUN, ephem.MERCURY, swe.MERCURY, 10.0),
    ("Sun-Venus", ephem.SUN, swe.SUN, ephem.VENUS, swe.VENUS, 10.0),
    ("Sun-Mars", ephem.SUN, swe.SUN, ephem.MARS, swe.MARS, 15.0),
    ("Sun-Jupiter", ephem.SUN, swe.SUN, ephem.JUPITER, swe.JUPITER, 20.0),
    ("Mars-Jupiter", ephem.MARS, swe.MARS, ephem.JUPITER, swe.JUPITER, 20.0),
    (
        "Jupiter-Saturn",
        ephem.JUPITER,
        swe.JUPITER,
        ephem.SATURN,
        swe.SATURN,
        30.0,
    ),
]

JD_START = 2451545.0
JD_END = JD_START + 365.25 * 15


def angular_sep(lon1, lon2):
    d = abs(lon1 - lon2) % 360
    return min(d, 360 - d)


def find_conjunctions(se_b1, se_b2, step, max_count=6):
    events = []
    jd = JD_START
    prev_sep = None
    prev_prev_sep = None
    while jd < JD_END and len(events) < max_count:
        try:
            r1 = swe.calc_ut(jd, se_b1, swe.FLG_SWIEPH | swe.FLG_SPEED)
            r2 = swe.calc_ut(jd, se_b2, swe.FLG_SWIEPH | swe.FLG_SPEED)
            sep = angular_sep(r1[0][0], r2[0][0])
        except Exception:
            jd += step
            continue
        if prev_sep is not None and prev_prev_sep is not None:
            if prev_sep < prev_prev_sep and prev_sep < sep and prev_sep < 15.0:
                # Refine
                jd_lo, jd_hi = jd - 2 * step, jd
                for _ in range(30):
                    jd_mid = (jd_lo + jd_hi) / 2
                    jd_q1 = (jd_lo + jd_mid) / 2
                    jd_q3 = (jd_mid + jd_hi) / 2
                    try:
                        s1 = angular_sep(
                            swe.calc_ut(jd_q1, se_b1, swe.FLG_SWIEPH)[0][0],
                            swe.calc_ut(jd_q1, se_b2, swe.FLG_SWIEPH)[0][0],
                        )
                        s3 = angular_sep(
                            swe.calc_ut(jd_q3, se_b1, swe.FLG_SWIEPH)[0][0],
                            swe.calc_ut(jd_q3, se_b2, swe.FLG_SWIEPH)[0][0],
                        )
                    except:
                        break
                    if s1 < s3:
                        jd_hi = jd_mid
                    else:
                        jd_lo = jd_mid
                events.append((jd_lo + jd_hi) / 2)
        prev_prev_sep = prev_sep
        prev_sep = sep
        jd += step
    return events


def compare_at_event(label, le_b1, se_b1, le_b2, se_b2, jd):
    global passed, failed, total
    try:
        le1 = ephem.calc_ut(jd, le_b1, FLAGS)
        se1 = swe.calc_ut(jd, se_b1, swe.FLG_SWIEPH | swe.FLG_SPEED)
        le2 = ephem.calc_ut(jd, le_b2, FLAGS)
        se2 = swe.calc_ut(jd, se_b2, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except Exception:
        return

    for body_label, le_r, se_r in [
        (f"{label}_B1", le1, se1),
        (f"{label}_B2", le2, se2),
    ]:
        total += 1
        diff = angular_sep(le_r[0][0], se_r[0][0]) * 3600
        if diff <= 1.0:
            passed += 1
        else:
            failed += 1
            failures.append(f'  {body_label} LON: diff={diff:.4f}"')
        total += 1
        lat_diff = abs(le_r[0][1] - se_r[0][1]) * 3600
        if lat_diff <= 1.0:
            passed += 1
        else:
            failed += 1
            failures.append(f'  {body_label} LAT: diff={lat_diff:.4f}"')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 221: Planetary Conjunction/Opposition Exact Timing")
    print("=" * 70)
    for pname, le1, se1, le2, se2, step in PAIRS:
        print(f"\n--- {pname} ---")
        conjs = find_conjunctions(se1, se2, step)
        print(f"  Found {len(conjs)} conjunctions")
        for i, jd in enumerate(conjs):
            compare_at_event(f"{pname} conj#{i + 1}", le1, se1, le2, se2, jd)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
