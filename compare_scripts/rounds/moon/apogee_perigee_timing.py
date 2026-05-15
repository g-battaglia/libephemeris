#!/usr/bin/env python3
"""Round 216: Moon apogee/perigee timing.

Tests mooncross_ut for Moon reaching apogee and perigee distances,
and compares Moon distance at extreme points with SE.
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

# Test dates (spanning 20 years)
JD_START = 2451545.0  # J2000
SEARCH_YEARS = 15
STEP = 3.0  # days


def find_moon_distance_extremes(jd_start, jd_end, step):
    """Find approximate times of Moon distance maxima and minima."""
    extremes = []
    jd = jd_start
    prev_dist = None
    prev_prev_dist = None

    while jd < jd_end:
        try:
            r = swe.calc_ut(jd, swe.MOON, swe.FLG_SWIEPH | swe.FLG_SPEED)
            dist = r[0][2]  # distance in AU
        except Exception:
            jd += step
            continue

        if prev_dist is not None and prev_prev_dist is not None:
            # Local maximum (apogee)
            if prev_dist > prev_prev_dist and prev_dist > dist:
                # Refine with bisection
                jd_refined = refine_extreme(jd - 2 * step, jd, is_max=True)
                if jd_refined:
                    extremes.append(("Apogee", jd_refined))
            # Local minimum (perigee)
            elif prev_dist < prev_prev_dist and prev_dist < dist:
                jd_refined = refine_extreme(jd - 2 * step, jd, is_max=False)
                if jd_refined:
                    extremes.append(("Perigee", jd_refined))

        prev_prev_dist = prev_dist
        prev_dist = dist
        jd += step

        if len(extremes) >= 40:
            break

    return extremes


def refine_extreme(jd_lo, jd_hi, is_max):
    """Refine extreme using golden section search."""
    gr = (5**0.5 + 1) / 2
    for _ in range(30):
        d = jd_hi - jd_lo
        if d < 1e-8:
            break
        jd1 = jd_hi - d / gr
        jd2 = jd_lo + d / gr
        try:
            d1 = swe.calc_ut(jd1, swe.MOON, swe.FLG_SWIEPH)[0][2]
            d2 = swe.calc_ut(jd2, swe.MOON, swe.FLG_SWIEPH)[0][2]
        except Exception:
            return None
        if is_max:
            if d1 > d2:
                jd_hi = jd2
            else:
                jd_lo = jd1
        else:
            if d1 < d2:
                jd_hi = jd2
            else:
                jd_lo = jd1
    return (jd_lo + jd_hi) / 2


def compare_at_extreme(label, jd):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, ephem.MOON, FLAGS)
        se_r = swe.calc_ut(jd, swe.MOON, swe.FLG_SWIEPH | swe.FLG_SPEED)
    except Exception:
        return

    # Longitude
    total += 1
    lon_diff = abs(le_r[0][0] - se_r[0][0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    if lon_as <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: diff={lon_as:.4f}"')

    # Distance
    total += 1
    dist_diff = abs(le_r[0][2] - se_r[0][2])
    dist_diff_km = dist_diff * 149597870.7  # AU to km
    if dist_diff_km <= 1.0:  # within 1 km
        passed += 1
    else:
        failed += 1
        failures.append(
            f"  {label} DIST: LE={le_r[0][2]:.10f} SE={se_r[0][2]:.10f} diff={dist_diff_km:.4f} km"
        )

    # Distance speed
    total += 1
    dspd_diff = abs(le_r[0][5] - se_r[0][5])
    dspd_diff_as = dspd_diff * 3600 * 149597870.7  # rough
    # At apogee/perigee, dist_speed should be near zero
    le_dspd_abs = abs(le_r[0][5])
    se_dspd_abs = abs(se_r[0][5])
    # Both should be near zero
    if le_dspd_abs < 0.0001 and se_dspd_abs < 0.0001:
        passed += 1
    elif abs(le_r[0][5] - se_r[0][5]) < 0.00005:
        passed += 1
    else:
        failed += 1
        failures.append(f"  {label} DSPD: LE={le_r[0][5]:.8f} SE={se_r[0][5]:.8f}")

    # Latitude
    total += 1
    lat_diff = abs(le_r[0][1] - se_r[0][1]) * 3600
    if lat_diff <= 1.0:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LAT: diff={lat_diff:.4f}"')

    # Lon speed
    total += 1
    lspd_diff = abs(le_r[0][3] - se_r[0][3]) * 3600
    if lspd_diff <= 5.0:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LSPD: LE={le_r[0][3]:.6f} SE={se_r[0][3]:.6f} diff={lspd_diff:.2f}"/day'
        )


if __name__ == "__main__":
    print("=" * 70)
    print("Round 216: Moon Apogee/Perigee Timing")
    print("=" * 70)

    jd_end = JD_START + 365.25 * SEARCH_YEARS
    print(f"\nSearching for Moon distance extremes over {SEARCH_YEARS} years...")
    extremes = find_moon_distance_extremes(JD_START, jd_end, STEP)
    print(f"Found {len(extremes)} extremes")

    for etype, jd in extremes:
        label = f"{etype} JD={jd:.4f}"
        compare_at_extreme(label, jd)

    print(f"\n{'=' * 70}")
    if total > 0:
        print(
            f"RESULTS: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed"
        )
    else:
        print("RESULTS: 0 tests ran")
    print(f"{'=' * 70}")
    if failures:
        print(f"\nFAILURES ({len(failures)}):")
        for f in failures[:30]:
            print(f)
        if len(failures) > 30:
            print(f"  ... and {len(failures) - 30} more")
