#!/usr/bin/env python3
"""Round 225: Equinox/solstice Sun positions.

Tests Sun longitude at equinox (0°/180°) and solstice (90°/270°) times,
plus nutation, obliquity, and sidereal time at those critical moments.
"""

from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED
JD_START = 2451545.0


def find_equinoxes_solstices(years=25):
    """Find approximate equinox/solstice times."""
    events = []
    jd = JD_START
    step = 1.0
    prev_lon = None
    while jd < JD_START + 365.25 * years:
        try:
            r = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH)
            lon = r[0][0]
        except:
            jd += step
            continue

        if prev_lon is not None:
            for target in [0.0, 90.0, 180.0, 270.0]:
                crossed = False
                if target == 0.0:
                    crossed = prev_lon > 350 and lon < 10
                else:
                    crossed = prev_lon < target and lon >= target
                if crossed:
                    # Refine
                    jd_lo, jd_hi = jd - step, jd
                    for _ in range(30):
                        mid = (jd_lo + jd_hi) / 2
                        try:
                            ml = swe.calc_ut(mid, swe.SUN, swe.FLG_SWIEPH)[0][0]
                        except:
                            break
                        if target == 0.0:
                            if ml > 180:
                                jd_lo = mid
                            else:
                                jd_hi = mid
                        else:
                            if ml < target:
                                jd_lo = mid
                            else:
                                jd_hi = mid
                    etype = {
                        0: "VernalEq",
                        90: "SummerSol",
                        180: "AutumnalEq",
                        270: "WinterSol",
                    }[target]
                    events.append((etype, (jd_lo + jd_hi) / 2))

        prev_lon = lon
        jd += step
        if len(events) >= 100:
            break
    return events


def compare_at(label, jd):
    global passed, failed, total

    # Sun position
    try:
        le_s = ephem.calc_ut(jd, ephem.SUN, FLAGS)
        se_s = swe.calc_ut(jd, swe.SUN, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except:
        return

    total += 1
    d = abs(le_s[0][0] - se_s[0][0])
    if d > 180:
        d = 360 - d
    if d * 3600 <= 0.5:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} Sun LON: diff={d * 3600:.4f}"')

    # Obliquity
    try:
        le_n = ephem.calc_ut(jd, -1, 0)
        se_n = swe.calc_ut(jd, -1, swe.FLG_SWIEPH)
    except:
        return

    total += 1
    d = abs(le_n[0][0] - se_n[0][0]) * 3600
    if d <= 0.1:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} OBL: diff={d:.6f}"')

    # Nutation dpsi
    total += 1
    d = abs(le_n[0][2] - se_n[0][2]) * 3600
    if d <= 0.1:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} DPSI: diff={d:.6f}"')

    # Sidereal time
    try:
        le_st = ephem.sidtime(jd)
        se_st = swe.sidtime(jd)
    except:
        return

    total += 1
    d = abs(le_st - se_st)
    if d > 12:
        d = 24 - d
    if d * 3600 <= 0.1:
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label} SIDTIME: diff={d * 3600:.4f}s")


if __name__ == "__main__":
    print("=" * 70)
    print("Round 225: Equinox/Solstice Sun Positions")
    print("=" * 70)
    events = find_equinoxes_solstices(25)
    print(f"Found {len(events)} equinoxes/solstices")
    for etype, jd in events:
        compare_at(f"{etype} JD={jd:.4f}", jd)

    print(f"\n{'=' * 70}")
    if total > 0:
        print(
            f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
        )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
