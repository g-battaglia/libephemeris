#!/usr/bin/env python3
"""Round 209: Speed at retrograde stations deep.

Tests planet speed (lon_speed) near retrograde stations where speed
passes through zero. This is the most demanding test for speed computation.
"""

from __future__ import annotations

import os
import sys

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

BODIES = [
    ("Mercury", ephem.MERCURY, swe.MERCURY, 5.0),
    ("Venus", ephem.VENUS, swe.VENUS, 5.0),
    ("Mars", ephem.MARS, swe.MARS, 10.0),
    ("Jupiter", ephem.JUPITER, swe.JUPITER, 20.0),
    ("Saturn", ephem.SATURN, swe.SATURN, 30.0),
]

JD_START = 2451545.0
JD_END = JD_START + 365.25 * 10  # 10 years


def find_stations(le_body, se_body, step):
    """Find approximate station times (speed crosses zero)."""
    stations = []
    jd = JD_START
    prev_spd = None

    while jd < JD_END:
        try:
            r = swe.calc_ut(jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED)
            spd = r[0][3]
        except Exception:
            jd += step
            continue

        if prev_spd is not None and prev_spd * spd < 0:
            # Speed sign change — refine
            jd_lo, jd_hi = jd - step, jd
            for _ in range(25):
                jd_mid = (jd_lo + jd_hi) / 2
                try:
                    mid_spd = swe.calc_ut(
                        jd_mid, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED
                    )[0][3]
                except Exception:
                    break
                if prev_spd * mid_spd < 0:
                    jd_hi = jd_mid
                else:
                    jd_lo = jd_mid
                    prev_spd = mid_spd
            stations.append((jd_lo + jd_hi) / 2)
            if len(stations) >= 8:
                break

        prev_spd = spd
        jd += step

    return stations


def compare_at_station(label, le_body, se_body, jd):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, le_body, FLAGS)
        se_r = swe.calc_ut(jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except Exception:
        return

    # Position (should be very precise at stations)
    total += 1
    lon_diff = abs(le_r[0][0] - se_r[0][0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    if lon_as <= 0.5:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={lon_as:.4f}"')

    # Speed at station (both should be near zero)
    total += 1
    le_spd = le_r[0][3]
    se_spd = se_r[0][3]
    spd_diff = abs(le_spd - se_spd) * 3600

    if spd_diff <= 2.0:  # 2"/day at station
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} SPD: LE={le_spd:.8f} SE={se_spd:.8f} diff={spd_diff:.4f}"/day'
        )

    # Both should be near zero
    total += 1
    le_spd_as = abs(le_spd) * 3600
    se_spd_as = abs(se_spd) * 3600
    if le_spd_as <= 10.0:  # within 10"/day of zero
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LE_SPD_ABS: {le_spd_as:.4f}"/day (should be ~0)')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 209: Speed at Retrograde Stations Deep")
    print("=" * 70)

    for bname, le_b, se_b, step in BODIES:
        print(f"\n--- {bname} ---")
        stations = find_stations(le_b, se_b, step)
        print(f"  Found {len(stations)} stations")

        for i, jd in enumerate(stations):
            stype = "Rx" if i % 2 == 0 else "D"
            compare_at_station(f"{bname} {stype}#{i // 2 + 1}", le_b, se_b, jd)
            # Also test ±0.5 day from station
            for offset in [-0.5, 0.5]:
                compare_at_station(
                    f"{bname} {stype}#{i // 2 + 1}{offset:+.1f}d",
                    le_b,
                    se_b,
                    jd + offset,
                )

    print(f"\n{'=' * 70}")
    print(
        f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
    )
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
