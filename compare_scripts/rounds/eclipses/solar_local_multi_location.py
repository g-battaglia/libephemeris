#!/usr/bin/env python3
"""Round 48: Solar Eclipse Local (sol_eclipse_when_loc) Multi-Location Sweep.

Tests solar eclipse local predictions across multiple geographic locations
and eclipse events, comparing libephemeris against pyswisseph.

Phases:
  P1: Known total/annular eclipses at path-centerline locations
  P2: Partial eclipses at off-path locations
  P3: sol_eclipse_how at known eclipse maximum times
  P4: Eclipse attributes (magnitude, obscuration, diameter ratio)
  P5: Sunrise/sunset eclipses (horizon events)
  P6: Southern hemisphere and equatorial eclipses
"""

from __future__ import annotations

import math
import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

passed = 0
failed = 0
errors = 0
results = {"passed": [], "failed": [], "errors": []}

FLG_SWIEPH = 2


def record(phase, label, ok, detail=""):
    global passed, failed
    if ok:
        passed += 1
        results["passed"].append(f"{phase} {label}: {detail}")
    else:
        failed += 1
        results["failed"].append(f"{phase} {label}: {detail}")


def safe_eclipse_when_loc(jd_start, lon, lat, alt=0.0, backward=False):
    """Get sol_eclipse_when_loc from both SE and LE."""
    geopos = [lon, lat, alt]
    try:
        # pyswisseph: sol_eclipse_when_loc(tjdut, geopos, flags, backwards)
        se_ret = swe.sol_eclipse_when_loc(jd_start, geopos, FLG_SWIEPH, backward)
        se_flag = se_ret[0]
        se_tret = se_ret[1]
        se_attr = se_ret[2]
    except Exception as e:
        return None, f"SE error: {e}"

    try:
        # libephemeris: sol_eclipse_when_loc(tjd_start, ifl, geopos, backward)
        le_ret = ephem.sol_eclipse_when_loc(
            jd_start, FLG_SWIEPH, geopos, backward
        )
        le_flag = le_ret[0]
        le_tret = le_ret[1]
        le_attr = le_ret[2]
    except Exception as e:
        return None, f"LE error: {e}"

    return (se_flag, se_tret, se_attr, le_flag, le_tret, le_attr), ""


def safe_eclipse_how(jd, lon, lat, alt=0.0):
    """Get sol_eclipse_how from both SE and LE."""
    geopos = [lon, lat, alt]
    try:
        # pyswisseph: sol_eclipse_how(tjdut, geopos, flags)
        se_ret = swe.sol_eclipse_how(jd, geopos, FLG_SWIEPH)
        se_flag = se_ret[0]
        se_attr = se_ret[1]
    except Exception as e:
        return None, f"SE error: {e}"

    try:
        # libephemeris: sol_eclipse_how(tjd_ut, ifl, geopos)
        le_ret = ephem.sol_eclipse_how(jd, FLG_SWIEPH, geopos)
        le_flag = le_ret[0]
        le_attr = le_ret[1]
    except Exception as e:
        return None, f"LE error: {e}"

    return (se_flag, se_attr, le_flag, le_attr), ""


def compare_eclipse_times(
    phase, label, se_tret, le_tret, max_tol_s=120, contact_tol_s=300
):
    """Compare eclipse times. Returns True if all within tolerance."""
    global errors
    time_names = ["max", "C1", "C2", "C3", "C4", "sunrise", "sunset"]
    worst_name = ""
    worst_diff = 0.0
    all_ok = True
    details = []

    for i, name in enumerate(time_names):
        if i >= len(se_tret) or i >= len(le_tret):
            break
        se_t = se_tret[i]
        le_t = le_tret[i]

        # Skip if both are 0 (event doesn't occur)
        if se_t == 0.0 and le_t == 0.0:
            continue

        # If one is 0 and other isn't, that's a structural difference
        if (se_t == 0.0) != (le_t == 0.0):
            details.append(f"{name}:SE={se_t:.6f}/LE={le_t:.6f}(MISSING)")
            # Only fail for max, C1, C4 (critical). C2/C3 can differ for near-limit events.
            if name in ("max", "C1", "C4"):
                all_ok = False
                worst_diff = 999999
                worst_name = name
            continue

        diff_s = abs(se_t - le_t) * 86400  # convert days to seconds
        tol = max_tol_s if name == "max" else contact_tol_s

        if diff_s > tol:
            all_ok = False

        if diff_s > worst_diff:
            worst_diff = diff_s
            worst_name = name

        if diff_s > 10:  # Only show non-trivial diffs
            details.append(f"{name}:{diff_s:.1f}s")

    detail_str = f"worst={worst_name}:{worst_diff:.1f}s " + " ".join(details)
    record(phase, label, all_ok, detail_str)
    return all_ok


def compare_eclipse_attrs(
    phase,
    label,
    se_attr,
    le_attr,
    mag_tol=0.02,
    obs_tol=0.03,
    ratio_tol=0.01,
    alt_tol=1.0,
    az_tol=2.0,
    width_tol=50.0,
):
    """Compare eclipse attributes."""
    all_ok = True
    details = []

    # attr[0]: magnitude
    diff_mag = abs(se_attr[0] - le_attr[0])
    if diff_mag > mag_tol:
        all_ok = False
    if diff_mag > 0.001:
        details.append(f"mag:SE={se_attr[0]:.4f}/LE={le_attr[0]:.4f}/d={diff_mag:.4f}")

    # attr[1]: diameter ratio
    diff_ratio = abs(se_attr[1] - le_attr[1])
    if diff_ratio > ratio_tol:
        all_ok = False
    if diff_ratio > 0.001:
        details.append(
            f"ratio:SE={se_attr[1]:.4f}/LE={le_attr[1]:.4f}/d={diff_ratio:.4f}"
        )

    # attr[2]: obscuration
    diff_obs = abs(se_attr[2] - le_attr[2])
    if diff_obs > obs_tol:
        all_ok = False
    if diff_obs > 0.001:
        details.append(f"obs:SE={se_attr[2]:.4f}/LE={le_attr[2]:.4f}/d={diff_obs:.4f}")

    # attr[3]: shadow width (km)
    if se_attr[3] != 0 or le_attr[3] != 0:
        diff_w = abs(se_attr[3] - le_attr[3])
        rel_w = diff_w / max(abs(se_attr[3]), 1.0)
        if diff_w > width_tol and rel_w > 0.25:
            all_ok = False
        if diff_w > 1.0:
            details.append(
                f"width:SE={se_attr[3]:.1f}/LE={le_attr[3]:.1f}/d={diff_w:.1f}km"
            )

    # attr[4]: azimuth (skip if both ~0)
    if abs(se_attr[4]) > 0.1 or abs(le_attr[4]) > 0.1:
        diff_az = abs(se_attr[4] - le_attr[4])
        if diff_az > 180:
            diff_az = 360 - diff_az
        if diff_az > az_tol:
            all_ok = False
        if diff_az > 0.1:
            details.append(f"az:d={diff_az:.2f}°")

    # attr[5]: true altitude
    if abs(se_attr[5]) > 0.01 or abs(le_attr[5]) > 0.01:
        diff_alt = abs(se_attr[5] - le_attr[5])
        if diff_alt > alt_tol:
            all_ok = False
        if diff_alt > 0.05:
            details.append(
                f"alt:SE={se_attr[5]:.2f}/LE={le_attr[5]:.2f}/d={diff_alt:.2f}°"
            )

    detail_str = " ".join(details) if details else "all match"
    record(phase, label, all_ok, detail_str)
    return all_ok


# ── Eclipse events database ──
# Known solar eclipses with approximate centerline locations
ECLIPSES = {
    # Total solar eclipses
    "Total_2024Apr08_Texas": {
        "jd_before": 2460380.0,  # ~10 days before
        "lon": -98.5,
        "lat": 31.8,
        "type": "total",
    },
    "Total_2024Apr08_Maine": {
        "jd_before": 2460380.0,
        "lon": -69.0,
        "lat": 46.5,
        "type": "total",
    },
    "Total_2017Aug21_Oregon": {
        "jd_before": 2457975.0,
        "lon": -121.2,
        "lat": 44.6,
        "type": "total",
    },
    "Total_2017Aug21_SCarolina": {
        "jd_before": 2457975.0,
        "lon": -79.9,
        "lat": 34.0,
        "type": "total",
    },
    "Total_2019Jul02_Chile": {
        "jd_before": 2458660.0,
        "lon": -70.5,
        "lat": -30.5,
        "type": "total",
    },
    # Annular solar eclipses
    "Annular_2023Oct14_Texas": {
        "jd_before": 2460215.0,
        "lon": -100.4,
        "lat": 29.5,
        "type": "annular",
    },
    "Annular_2023Oct14_Yucatan": {
        "jd_before": 2460215.0,
        "lon": -88.5,
        "lat": 20.5,
        "type": "annular",
    },
    "Annular_2020Jun21_India": {
        "jd_before": 2459010.0,
        "lon": 79.0,
        "lat": 23.0,
        "type": "annular",
    },
    # Partial eclipses (off center)
    "Partial_2024Apr08_Chicago": {
        "jd_before": 2460380.0,
        "lon": -87.6,
        "lat": 41.9,
        "type": "partial",
    },
    "Partial_2024Apr08_NYC": {
        "jd_before": 2460380.0,
        "lon": -74.0,
        "lat": 40.7,
        "type": "partial",
    },
    "Partial_2024Apr08_London": {
        "jd_before": 2460380.0,
        "lon": -0.1,
        "lat": 51.5,
        "type": "partial",
    },
    # Eclipses at various latitudes
    "Total_2021Dec04_Antarctica": {
        "jd_before": 2459545.0,
        "lon": 0.0,
        "lat": -75.0,
        "type": "total",
    },
    "Annular_2021Jun10_Canada": {
        "jd_before": 2459368.0,
        "lon": -80.0,
        "lat": 50.0,
        "type": "annular",
    },
}


def phase1():
    """Known eclipses at/near centerline — timing comparison."""
    global errors
    print("\n=== P1: Known eclipses timing ===")

    for name, info in ECLIPSES.items():
        try:
            result, err = safe_eclipse_when_loc(
                info["jd_before"], info["lon"], info["lat"]
            )
            if result is None:
                errors += 1
                results["errors"].append(f"P1 {name}: {err}")
                continue

            se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

            compare_eclipse_times(
                "P1", name, se_tret, le_tret, max_tol_s=120, contact_tol_s=300
            )
        except Exception as e:
            errors += 1
            results["errors"].append(f"P1 {name}: {e}")


def phase2():
    """Partial eclipses at off-path locations."""
    global errors
    print("\n=== P2: Partial eclipse off-path timing ===")

    # Cities that see partial eclipse during 2024 Apr 08
    partial_locs = {
        "Miami": (-80.2, 25.8),
        "Denver": (-104.9, 39.7),
        "Toronto": (-79.4, 43.7),
        "Montreal": (-73.6, 45.5),
        "Washington": (-77.0, 38.9),
        "Atlanta": (-84.4, 33.7),
        "Havana": (-82.4, 23.1),
    }

    jd_before = 2460380.0
    for city_name, (lon, lat) in partial_locs.items():
        try:
            result, err = safe_eclipse_when_loc(jd_before, lon, lat)
            if result is None:
                errors += 1
                results["errors"].append(f"P2 {city_name}: {err}")
                continue

            se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

            # Check that both found the same eclipse
            if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                record("P2", f"{city_name} no_eclipse", True, "both: no eclipse")
                continue

            compare_eclipse_times(
                "P2",
                f"{city_name}_times",
                se_tret,
                le_tret,
                max_tol_s=120,
                contact_tol_s=300,
            )
        except Exception as e:
            errors += 1
            results["errors"].append(f"P2 {city_name}: {e}")


def phase3():
    """sol_eclipse_how at known eclipse maximum times."""
    global errors
    print("\n=== P3: sol_eclipse_how at eclipse maximum ===")

    # First, get eclipse max times from SE, then test how() at those times
    test_cases = {
        "2024Apr08_Texas": (2460380.0, -98.5, 31.8),
        "2024Apr08_NYC": (2460380.0, -74.0, 40.7),
        "2017Aug21_Oregon": (2457975.0, -121.2, 44.6),
        "2023Oct14_Texas": (2460215.0, -100.4, 29.5),
        "2019Jul02_Chile": (2458660.0, -70.5, -30.5),
    }

    for name, (jd_before, lon, lat) in test_cases.items():
        try:
            # Get eclipse max time from SE
            geopos = [lon, lat, 0.0]
            se_when = swe.sol_eclipse_when_loc(jd_before, geopos, FLG_SWIEPH, False)
            jd_max = se_when[1][0]  # tret[0] = maximum

            if jd_max == 0.0:
                record("P3", f"{name}_how", True, "no eclipse found by SE")
                continue

            # Now test how() at that time
            result, err = safe_eclipse_how(jd_max, lon, lat)
            if result is None:
                errors += 1
                results["errors"].append(f"P3 {name}: {err}")
                continue

            se_flag, se_attr, le_flag, le_attr = result

            # Compare attributes
            compare_eclipse_attrs("P3", f"{name}_how", se_attr, le_attr)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P3 {name}: {e}")


def phase4():
    """Eclipse attributes comparison (magnitude, obscuration, ratio)."""
    global errors
    print("\n=== P4: Eclipse attributes ===")

    for name, info in ECLIPSES.items():
        try:
            result, err = safe_eclipse_when_loc(
                info["jd_before"], info["lon"], info["lat"]
            )
            if result is None:
                continue

            se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

            if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                continue

            compare_eclipse_attrs("P4", f"{name}_attr", se_attr, le_attr)
        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 {name}: {e}")


def phase5():
    """Sunrise/sunset eclipses — locations near eclipse path edges."""
    global errors
    print("\n=== P5: Horizon eclipses ===")

    # Locations where eclipse may be near sunrise/sunset
    horizon_cases = {
        # 2024 Apr 08 — sunrise eclipse in western US
        "2024Apr08_SanFran": (2460380.0, -122.4, 37.8),
        # 2017 Aug 21 — sunset eclipse in eastern Atlantic
        "2017Aug21_Bermuda": (2457975.0, -64.8, 32.3),
        # 2023 Oct 14 — partial at sunset for eastern Europe
        "2023Oct14_Portugal": (2460215.0, -9.1, 38.7),
        # 2020 Jun 21 annular — sunrise in Africa
        "2020Jun21_Congo": (2459010.0, 25.0, 0.0),
    }

    for name, (jd_before, lon, lat) in horizon_cases.items():
        try:
            result, err = safe_eclipse_when_loc(jd_before, lon, lat)
            if result is None:
                errors += 1
                results["errors"].append(f"P5 {name}: {err}")
                continue

            se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

            if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                record("P5", f"{name} no_eclipse", True, "both: no eclipse visible")
                continue

            # Wider tolerance for horizon events (refraction modeling differs)
            compare_eclipse_times(
                "P5",
                f"{name}_times",
                se_tret,
                le_tret,
                max_tol_s=180,
                contact_tol_s=600,
            )
        except Exception as e:
            errors += 1
            results["errors"].append(f"P5 {name}: {e}")


def phase6():
    """Southern hemisphere and equatorial eclipses."""
    global errors
    print("\n=== P6: Southern/equatorial eclipses ===")

    south_cases = {
        # 2024 Oct 02 annular — southern Pacific/Chile/Argentina
        "2024Oct02_Santiago": (2460570.0, -70.7, -33.4),
        "2024Oct02_BuenosAires": (2460570.0, -58.4, -34.6),
        # 2028 Jul 22 total — Australia/NZ
        "2028Jul22_Sydney": (2461930.0, 151.2, -33.9),
        # 2026 Aug 12 total — Iceland to Iberia
        "2026Aug12_Madrid": (2461220.0, -3.7, 40.4),
        "2026Aug12_Reykjavik": (2461220.0, -22.0, 64.1),
        # 2030 Nov 25 total — southern Africa
        "2030Nov25_CapeTown": (2462770.0, 18.4, -33.9),
        "2030Nov25_Johannesburg": (2462770.0, 28.0, -26.2),
        # Equatorial
        "2019Dec26_Singapore": (2458843.0, 103.8, 1.3),
        "2019Dec26_Manila": (2458843.0, 121.0, 14.6),
    }

    for name, (jd_before, lon, lat) in south_cases.items():
        try:
            result, err = safe_eclipse_when_loc(jd_before, lon, lat)
            if result is None:
                errors += 1
                results["errors"].append(f"P6 {name}: {err}")
                continue

            se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

            if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                record("P6", f"{name} no_eclipse", True, "both: no eclipse found")
                continue

            if se_tret[0] == 0.0 or le_tret[0] == 0.0:
                record(
                    "P6",
                    f"{name}",
                    False,
                    f"one found eclipse, other didn't: SE_max={se_tret[0]:.2f} LE_max={le_tret[0]:.2f}",
                )
                continue

            # Check they found the same eclipse (within 30 days)
            if abs(se_tret[0] - le_tret[0]) > 30:
                record(
                    "P6",
                    f"{name}",
                    False,
                    f"different eclipses: SE={se_tret[0]:.2f} LE={le_tret[0]:.2f}",
                )
                continue

            compare_eclipse_times(
                "P6",
                f"{name}_times",
                se_tret,
                le_tret,
                max_tol_s=120,
                contact_tol_s=300,
            )
            compare_eclipse_attrs("P6", f"{name}_attr", se_attr, le_attr)
        except Exception as e:
            errors += 1
            results["errors"].append(f"P6 {name}: {e}")


def main():
    print("=" * 70)
    print("ROUND 48: Solar Eclipse Local Multi-Location Sweep")
    print("=" * 70)

    phase1()
    print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

    phase2()
    print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

    phase3()
    print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

    phase4()
    print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

    phase5()
    print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

    phase6()
    print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

    total = passed + failed + errors
    pct = 100 * passed / total if total else 0
    print("\n" + "=" * 70)
    print(f"ROUND 48 FINAL: {passed}/{total} passed ({pct:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print("=" * 70)

    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:40]:
            print(f)

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:10]:
            print(e)

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
