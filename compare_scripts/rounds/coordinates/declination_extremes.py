#!/usr/bin/env python3
"""Round 206: Planetary declination extremes.

Tests equatorial declination at extreme values (max/min declination)
for planets, particularly when planets are out-of-bounds (beyond ±23.44°).
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
total = 0
failures = []

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED | ephem.FLG_EQUATORIAL

BODIES = [
    ("Sun", ephem.SUN, swe.SUN, 0.5),
    ("Moon", ephem.MOON, swe.MOON, 1.0),
    ("Mercury", ephem.MERCURY, swe.MERCURY, 0.5),
    ("Venus", ephem.VENUS, swe.VENUS, 0.5),
    ("Mars", ephem.MARS, swe.MARS, 0.5),
    ("Jupiter", ephem.JUPITER, swe.JUPITER, 0.5),
    ("Saturn", ephem.SATURN, swe.SATURN, 0.5),
    ("Uranus", ephem.URANUS, swe.URANUS, 0.5),
    ("Neptune", ephem.NEPTUNE, swe.NEPTUNE, 0.5),
    ("Pluto", ephem.PLUTO, swe.PLUTO, 1.0),
]

# Scan for declination extremes across 25 years
JD_START = 2451545.0
JD_END = JD_START + 365.25 * 25


def find_declination_extremes(le_body, se_body, step=5.0):
    """Find dates of maximum and minimum declination."""
    extremes = []
    jd = JD_START
    prev_dec = None
    prev_prev_dec = None

    while jd < JD_END:
        try:
            se_r = swe.calc_ut(
                jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
            )
            dec = se_r[0][1]
        except Exception:
            jd += step
            continue

        if prev_dec is not None and prev_prev_dec is not None:
            # Check for local max/min
            if (prev_dec > dec and prev_dec > prev_prev_dec) or (
                prev_dec < dec and prev_dec < prev_prev_dec
            ):
                extremes.append(jd - step)
                if len(extremes) >= 6:
                    break

        prev_prev_dec = prev_dec
        prev_dec = dec
        jd += step

    return extremes


def compare_at(label, le_body, se_body, jd, tol):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, le_body, FLAGS)
        se_r = swe.calc_ut(
            jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL
        )
    except Exception:
        return

    # RA
    total += 1
    ra_diff = abs(le_r[0][0] - se_r[0][0])
    if ra_diff > 180:
        ra_diff = 360 - ra_diff
    ra_as = ra_diff * 3600
    if ra_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} RA: diff={ra_as:.4f}"')

    # Dec
    total += 1
    dec_as = abs(le_r[0][1] - se_r[0][1]) * 3600
    if dec_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} DEC: LE={le_r[0][1]:.6f} SE={se_r[0][1]:.6f} diff={dec_as:.4f}"'
        )

    # Dec speed
    total += 1
    dspd = abs(le_r[0][4] - se_r[0][4]) * 3600
    if dspd <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} DEC_SPD: diff={dspd:.4f}"/day')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 206: Planetary Declination Extremes")
    print("=" * 70)

    for bname, le_b, se_b, tol in BODIES:
        step = 1.0 if bname == "Moon" else 10.0
        print(f"\n--- {bname} ---")
        extremes = find_declination_extremes(le_b, se_b, step)
        print(f"  Found {len(extremes)} declination extremes")

        for i, jd in enumerate(extremes):
            compare_at(f"{bname} ext#{i + 1}", le_b, se_b, jd, tol)
            compare_at(f"{bname} ext#{i + 1}-1d", le_b, se_b, jd - 1, tol)
            compare_at(f"{bname} ext#{i + 1}+1d", le_b, se_b, jd + 1, tol)

    # Also test solstice points (Sun at max/min declination)
    print("\n--- Sun at Solstices ---")
    solstice_jds = [
        2451723.5,
        2451907.5,  # 2000 summer/winter
        2455371.5,
        2455555.5,  # 2010
        2459022.5,
        2459206.5,  # 2020
        2460473.5,
        2460657.5,  # 2024
    ]
    for jd in solstice_jds:
        compare_at(f"Sun solstice JD={jd:.1f}", ephem.SUN, swe.SUN, jd, 0.5)

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
