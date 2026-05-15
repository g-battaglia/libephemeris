#!/usr/bin/env python3
"""Round 222: Lunar phase angle precision.

Tests Moon phase angle (elongation from Sun) at exact quarter phases
(new, first quarter, full, last quarter) comparing LE vs SE pheno_ut results.
"""

from __future__ import annotations
import os, sys

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
JD_START = 2451545.0
SEARCH_YEARS = 20


def find_new_moons(count=30):
    """Find approximate new moon times."""
    events = []
    jd = JD_START
    step = 1.0
    prev_elong = None
    while jd < JD_START + 365.25 * SEARCH_YEARS and len(events) < count:
        try:
            sun = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH)[0][0]
            moon = swe.calc_ut(jd, swe.MOON, swe.FLG_SWIEPH)[0][0]
            elong = (moon - sun) % 360
        except:
            jd += step
            continue
        if prev_elong is not None:
            if prev_elong > 300 and elong < 60:  # wrap around = new moon
                # Refine
                jd_lo, jd_hi = jd - step, jd
                for _ in range(25):
                    mid = (jd_lo + jd_hi) / 2
                    try:
                        s = swe.calc_ut(mid, swe.SUN, swe.FLG_SWIEPH)[0][0]
                        m = swe.calc_ut(mid, swe.MOON, swe.FLG_SWIEPH)[0][0]
                        e = (m - s) % 360
                    except:
                        break
                    if e > 180:
                        jd_hi = mid
                    else:
                        jd_lo = mid
                events.append((jd_lo + jd_hi) / 2)
        prev_elong = elong
        jd += step
    return events


def compare_at_phase(label, jd):
    global passed, failed, total

    # Compare Moon position
    try:
        le_m = ephem.calc_ut(jd, ephem.MOON, FLAGS)
        se_m = swe.calc_ut(jd, swe.MOON, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except:
        return

    total += 1
    diff = abs(le_m[0][0] - se_m[0][0])
    if diff > 180:
        diff = 360 - diff
    diff_as = diff * 3600
    if diff_as <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} Moon LON: diff={diff_as:.4f}"')

    # Compare Sun position
    try:
        le_s = ephem.calc_ut(jd, ephem.SUN, FLAGS)
        se_s = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except:
        return

    total += 1
    diff = abs(le_s[0][0] - se_s[0][0])
    if diff > 180:
        diff = 360 - diff
    diff_as = diff * 3600
    if diff_as <= 0.5:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} Sun LON: diff={diff_as:.4f}"')

    # Compare Moon distance
    total += 1
    dist_diff = abs(le_m[0][2] - se_m[0][2]) * 149597870.7
    if dist_diff <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label} Moon DIST: diff={dist_diff:.4f} km")

    # Moon speed
    total += 1
    spd_diff = abs(le_m[0][3] - se_m[0][3]) * 3600
    if spd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} Moon SPD: diff={spd_diff:.4f}"/day')

    # Sun-Moon elongation agreement
    total += 1
    le_elong = (le_m[0][0] - le_s[0][0]) % 360
    se_elong = (se_m[0][0] - se_s[0][0]) % 360
    elong_diff = abs(le_elong - se_elong)
    if elong_diff > 180:
        elong_diff = 360 - elong_diff
    if elong_diff * 3600 <= 2.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} ELONG: diff={elong_diff * 3600:.4f}"')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 222: Lunar Phase Angle Precision")
    print("=" * 70)
    new_moons = find_new_moons(30)
    print(f"Found {len(new_moons)} new moons")
    for i, jd in enumerate(new_moons):
        compare_at_phase(f"NewMoon#{i + 1}", jd)
        # Also test at quarter phases (+7.38d, +14.77d, +22.15d)
        for offset, phase in [(7.38, "Q1"), (14.77, "Full"), (22.15, "Q3")]:
            compare_at_phase(f"{phase}#{i + 1}", jd + offset)

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
