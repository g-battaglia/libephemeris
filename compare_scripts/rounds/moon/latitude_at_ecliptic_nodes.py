#!/usr/bin/env python3
"""Round 186: Planetary latitude at ecliptic nodes (lat≈0 crossing precision).

Tests planet positions near their ecliptic node crossings where latitude
passes through zero. This is a demanding test because small ephemeris
differences become most visible when the latitude is near zero.

We find approximate node crossing times by scanning, then compare
LE vs SE positions at those times and nearby offsets.
"""

from __future__ import annotations

import os
import sys
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

BODIES = {
    "Moon": (ephem.MOON, swe.MOON),
    "Mercury": (ephem.MERCURY, swe.MERCURY),
    "Venus": (ephem.VENUS, swe.VENUS),
    "Mars": (ephem.MARS, swe.MARS),
    "Jupiter": (ephem.JUPITER, swe.JUPITER),
    "Saturn": (ephem.SATURN, swe.SATURN),
}

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

# Test dates spanning 2000-2025
BASE_JD = 2451545.0  # J2000.0
SCAN_STEP = {
    "Moon": 1.0,  # Moon crosses ecliptic every ~13.6 days
    "Mercury": 5.0,
    "Venus": 10.0,
    "Mars": 20.0,
    "Jupiter": 50.0,
    "Saturn": 100.0,
}

passed = 0
failed = 0
total = 0
failures = []


def find_node_crossings(body_name, le_body, se_body, jd_start, jd_end, step):
    """Find approximate times when latitude crosses zero."""
    crossings = []
    jd = jd_start
    prev_lat = None
    while jd < jd_end:
        try:
            se_result = swe.calc_ut(jd, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED)
            se_lat = se_result[0][1]
        except Exception:
            jd += step
            continue

        if prev_lat is not None and prev_lat * se_lat < 0:
            # Sign change — refine with bisection
            jd_lo, jd_hi = jd - step, jd
            for _ in range(30):  # ~1e-9 day precision
                jd_mid = (jd_lo + jd_hi) / 2
                try:
                    mid_lat = swe.calc_ut(
                        jd_mid, se_body, swe.FLG_SWIEPH | swe.FLG_SPEED
                    )[0][1]
                except Exception:
                    break
                if prev_lat * mid_lat < 0:
                    jd_hi = jd_mid
                else:
                    jd_lo = jd_mid
                    prev_lat = mid_lat
            crossings.append((jd_lo + jd_hi) / 2)
            if len(crossings) >= 8:
                break

        prev_lat = se_lat
        jd += step
    return crossings


def compare_at_jd(body_name, le_body, se_body, jd, flags_le, flags_se, label):
    global passed, failed, total

    try:
        le_result = ephem.calc_ut(jd, le_body, flags_le)
        le_lon, le_lat, le_dist = le_result[0][0], le_result[0][1], le_result[0][2]
        le_lon_spd, le_lat_spd = le_result[0][3], le_result[0][4]
    except Exception as e:
        return

    try:
        se_result = swe.calc_ut(jd, se_body, flags_se)
        se_lon, se_lat, se_dist = se_result[0][0], se_result[0][1], se_result[0][2]
        se_lon_spd, se_lat_spd = se_result[0][3], se_result[0][4]
    except Exception:
        return

    # Longitude comparison
    total += 1
    lon_diff = abs(le_lon - se_lon)
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_diff_arcsec = lon_diff * 3600

    # Tolerance: Moon 1", others 0.5"
    lon_tol = 1.0 if body_name == "Moon" else 0.5
    if lon_diff_arcsec <= lon_tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LON: LE={le_lon:.8f} SE={se_lon:.8f} diff={lon_diff_arcsec:.4f}"'
        )

    # Latitude comparison (should be near zero at node crossing)
    total += 1
    lat_diff = abs(le_lat - se_lat)
    lat_diff_arcsec = lat_diff * 3600

    lat_tol = 1.0 if body_name == "Moon" else 0.5
    if lat_diff_arcsec <= lat_tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT: LE={le_lat:.8f} SE={se_lat:.8f} diff={lat_diff_arcsec:.4f}"'
        )

    # Latitude speed comparison
    total += 1
    lat_spd_diff = abs(le_lat_spd - se_lat_spd) * 3600

    lat_spd_tol = 5.0 if body_name == "Moon" else 2.0
    if lat_spd_tol == 0 or lat_spd_diff <= lat_spd_tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LAT_SPD: LE={le_lat_spd:.8f} SE={se_lat_spd:.8f} diff={lat_spd_diff:.4f}"/day'
        )

    # Longitude speed comparison
    total += 1
    lon_spd_diff = abs(le_lon_spd - se_lon_spd) * 3600

    lon_spd_tol = 5.0 if body_name == "Moon" else 2.0
    if lon_spd_diff <= lon_spd_tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f'  {label} LON_SPD: LE={le_lon_spd:.8f} SE={se_lon_spd:.8f} diff={lon_spd_diff:.4f}"/day'
        )


def test_at_node_crossings():
    """Test positions at and near ecliptic node crossings."""
    print("=" * 70)
    print("Round 186: Planetary Latitude at Ecliptic Nodes")
    print("=" * 70)

    jd_start = BASE_JD
    jd_end = BASE_JD + 365.25 * 25  # 25 years

    for body_name, (le_body, se_body) in BODIES.items():
        step = SCAN_STEP[body_name]
        print(f"\n--- {body_name} (scan step={step}d) ---")

        crossings = find_node_crossings(
            body_name, le_body, se_body, jd_start, jd_end, step
        )
        print(f"  Found {len(crossings)} node crossings")

        for i, jd_cross in enumerate(crossings):
            # Test at crossing
            compare_at_jd(
                body_name,
                le_body,
                se_body,
                jd_cross,
                FLAGS,
                swe.FLG_SWIEPH | swe.FLG_SPEED,
                f"{body_name} crossing#{i + 1} JD={jd_cross:.6f}",
            )

            # Test slightly before and after crossing (±0.1 day)
            for offset in [-0.1, 0.1]:
                compare_at_jd(
                    body_name,
                    le_body,
                    se_body,
                    jd_cross + offset,
                    FLAGS,
                    swe.FLG_SWIEPH | swe.FLG_SPEED,
                    f"{body_name} crossing#{i + 1}+{offset:+.1f}d",
                )

        # Also test with FLG_J2000
        print(f"  Testing J2000 frame at node crossings...")
        flags_le_j2k = FLAGS | ephem.FLG_J2000
        flags_se_j2k = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_J2000
        for i, jd_cross in enumerate(crossings[:4]):
            compare_at_jd(
                body_name,
                le_body,
                se_body,
                jd_cross,
                flags_le_j2k,
                flags_se_j2k,
                f"{body_name} J2000 crossing#{i + 1}",
            )

        # Test with FLG_NONUT
        print(f"  Testing NONUT at node crossings...")
        flags_le_nn = FLAGS | ephem.FLG_NONUT
        flags_se_nn = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_NONUT
        for i, jd_cross in enumerate(crossings[:4]):
            compare_at_jd(
                body_name,
                le_body,
                se_body,
                jd_cross,
                flags_le_nn,
                flags_se_nn,
                f"{body_name} NONUT crossing#{i + 1}",
            )


def test_equatorial_at_nodes():
    """Test equatorial coordinates at node crossings."""
    print(f"\n--- Equatorial at Node Crossings ---")
    global passed, failed, total

    jd_start = BASE_JD
    jd_end = BASE_JD + 365.25 * 10

    for body_name in ["Moon", "Mars", "Jupiter"]:
        le_body, se_body = BODIES[body_name]
        step = SCAN_STEP[body_name]
        crossings = find_node_crossings(
            body_name, le_body, se_body, jd_start, jd_end, step
        )

        flags_le = FLAGS | ephem.FLG_EQUATORIAL
        flags_se = swe.FLG_SWIEPH | swe.FLG_SPEED | swe.FLG_EQUATORIAL

        for i, jd_cross in enumerate(crossings[:4]):
            try:
                le_result = ephem.calc_ut(jd_cross, le_body, flags_le)
                se_result = swe.calc_ut(jd_cross, se_body, flags_se)

                # RA comparison
                total += 1
                ra_diff = abs(le_result[0][0] - se_result[0][0])
                if ra_diff > 180:
                    ra_diff = 360 - ra_diff
                ra_diff_arcsec = ra_diff * 3600

                tol = 1.0 if body_name == "Moon" else 0.5
                if ra_diff_arcsec <= tol:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  {body_name} EQ RA crossing#{i + 1}: diff={ra_diff_arcsec:.4f}"'
                    )

                # Dec comparison
                total += 1
                dec_diff = abs(le_result[0][1] - se_result[0][1]) * 3600
                if dec_diff <= tol:
                    passed += 1
                else:
                    failed += 1
                    failures.append(
                        f'  {body_name} EQ DEC crossing#{i + 1}: diff={dec_diff:.4f}"'
                    )

            except Exception as e:
                pass


if __name__ == "__main__":
    test_at_node_crossings()
    test_equatorial_at_nodes()

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
