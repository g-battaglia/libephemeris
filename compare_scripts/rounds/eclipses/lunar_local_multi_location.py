#!/usr/bin/env python3
"""Round 49: Lunar Eclipse Local (lun_eclipse_when_loc) Multi-Location Sweep.

Tests lunar eclipse local predictions across multiple geographic locations,
comparing libephemeris against pyswisseph.

Phases:
  P1: Known total/partial lunar eclipses at multiple locations
  P2: lun_eclipse_how at known eclipse maximum times
  P3: Eclipse attributes (magnitude, penumbral mag, altitude, azimuth)
  P4: Moonrise/moonset eclipses (horizon events)
  P5: Southern hemisphere and equatorial locations
  P6: Penumbral-only eclipses
"""

from __future__ import annotations

import math
import os
import sys

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

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


def safe_lun_when_loc(jd_start, lat, lon, alt=0.0):
    """Get lun_eclipse_when_loc from both SE and LE."""
    geopos_se = [lon, lat, alt]  # pyswisseph: [lon, lat, alt]
    try:
        se_ret = swe.lun_eclipse_when_loc(jd_start, geopos_se, FLG_SWIEPH, False)
        se_flag = se_ret[0]
        se_tret = se_ret[1]
        se_attr = se_ret[2]
    except Exception as e:
        return None, f"SE error: {e}"

    try:
        # libephemeris: lun_eclipse_when_loc(jd, geopos, flags)
        le_ret = ephem.lun_eclipse_when_loc(jd_start, (lon, lat, alt), FLG_SWIEPH)
        le_flag = le_ret[0]
        le_tret = le_ret[1]
        le_attr = le_ret[2]
    except Exception as e:
        return None, f"LE error: {e}"

    return (se_flag, se_tret, se_attr, le_flag, le_tret, le_attr), ""


def safe_lun_how(jd, lat, lon, alt=0.0):
    """Get lun_eclipse_how from both SE and LE."""
    geopos = [lon, lat, alt]
    try:
        # pyswisseph: lun_eclipse_how(tjd, geopos, ifl)
        se_ret = swe.lun_eclipse_how(jd, geopos, FLG_SWIEPH)
        se_flag = se_ret[0]
        se_attr = se_ret[1]
    except Exception as e:
        return None, f"SE error: {e}"

    try:
        # libephemeris: lun_eclipse_how(tjd, ifl, geopos)
        le_ret = ephem.lun_eclipse_how(jd, FLG_SWIEPH, geopos)
        le_flag = le_ret[0]
        le_attr = le_ret[1]
    except Exception as e:
        return None, f"LE error: {e}"

    return (se_flag, se_attr, le_flag, le_attr), ""


def compare_lun_times(phase, label, se_tret, le_tret, max_tol_s=120, phase_tol_s=300):
    """Compare lunar eclipse times."""
    # tret layout: [0]=max, [1]=reserved, [2]=partial_begin, [3]=partial_end,
    #              [4]=total_begin, [5]=total_end, [6]=penum_begin, [7]=penum_end,
    #              [8]=moonrise, [9]=moonset
    time_names = [
        "max",
        "rsv",
        "part_beg",
        "part_end",
        "tot_beg",
        "tot_end",
        "pen_beg",
        "pen_end",
        "moonrise",
        "moonset",
    ]
    worst_name = ""
    worst_diff = 0.0
    all_ok = True
    details = []

    for i, name in enumerate(time_names):
        if i >= len(se_tret) or i >= len(le_tret):
            break
        if i == 1:  # reserved
            continue

        se_t = se_tret[i]
        le_t = le_tret[i]

        if se_t == 0.0 and le_t == 0.0:
            continue

        if (se_t == 0.0) != (le_t == 0.0):
            details.append(f"{name}:MISSING(SE={se_t:.2f}/LE={le_t:.2f})")
            if name in ("max", "pen_beg", "pen_end"):
                all_ok = False
                worst_diff = 999999
                worst_name = name
            continue

        diff_s = abs(se_t - le_t) * 86400
        tol = max_tol_s if name == "max" else phase_tol_s

        if diff_s > tol:
            all_ok = False

        if diff_s > worst_diff:
            worst_diff = diff_s
            worst_name = name

        if diff_s > 10:
            details.append(f"{name}:{diff_s:.1f}s")

    detail_str = f"worst={worst_name}:{worst_diff:.1f}s " + " ".join(details)
    record(phase, label, all_ok, detail_str)
    return all_ok


def compare_lun_attrs(
    phase, label, se_attr, le_attr, mag_tol=0.02, alt_tol=1.0, az_tol=2.0
):
    """Compare lunar eclipse attributes."""
    all_ok = True
    details = []

    # attr[0]: umbral magnitude
    diff_mag = abs(se_attr[0] - le_attr[0])
    if diff_mag > mag_tol:
        all_ok = False
    if diff_mag > 0.001:
        details.append(f"umag:SE={se_attr[0]:.4f}/LE={le_attr[0]:.4f}/d={diff_mag:.4f}")

    # attr[1]: penumbral magnitude
    diff_pmag = abs(se_attr[1] - le_attr[1])
    if diff_pmag > mag_tol:
        all_ok = False
    if diff_pmag > 0.001:
        details.append(
            f"pmag:SE={se_attr[1]:.4f}/LE={le_attr[1]:.4f}/d={diff_pmag:.4f}"
        )

    # attr[4]: azimuth
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


# ── Known lunar eclipses ──
LOCATIONS = {
    "London": (51.5, -0.1),
    "NewYork": (40.7, -74.0),
    "LosAngeles": (34.1, -118.2),
    "Tokyo": (35.7, 139.7),
    "Sydney": (-33.9, 151.2),
    "CapeTown": (-33.9, 18.4),
    "Delhi": (28.6, 77.2),
    "SaoPaulo": (-23.5, -46.6),
    "Equator0": (0.0, 0.0),
    "Reykjavik": (64.1, -22.0),
}

LUNAR_ECLIPSES = {
    # Total lunar eclipses
    "Total_2022May16": 2459710.0,  # Total lunar eclipse
    "Total_2022Nov08": 2459890.0,  # Total lunar eclipse
    "Total_2025Mar14": 2460745.0,  # Total lunar eclipse
    "Total_2025Sep07": 2460923.0,  # Total lunar eclipse
    "Total_2018Jan31": 2458150.0,  # Total lunar eclipse
    "Total_2019Jan21": 2458505.0,  # Total lunar eclipse
    # Partial lunar eclipses
    "Partial_2024Sep18": 2460572.0,  # Partial
    "Partial_2023Oct28": 2460245.0,  # Partial
    "Partial_2021Nov19": 2459535.0,  # Partial
    # Penumbral
    "Penumbral_2024Mar25": 2460394.0,  # Penumbral
    "Penumbral_2023May05": 2460069.0,  # Penumbral
}


def phase1():
    """Known lunar eclipses at multiple locations — timing."""
    global errors
    print("\n=== P1: Known lunar eclipses timing ===")

    # Use a subset of locations for each eclipse
    test_locs = {
        k: v
        for k, v in LOCATIONS.items()
        if k in ("London", "NewYork", "Tokyo", "Sydney", "SaoPaulo")
    }

    for ecl_name, jd_before in LUNAR_ECLIPSES.items():
        for loc_name, (lat, lon) in test_locs.items():
            try:
                result, err = safe_lun_when_loc(jd_before, lat, lon)
                if result is None:
                    errors += 1
                    results["errors"].append(f"P1 {ecl_name} {loc_name}: {err}")
                    continue

                se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

                # Check both found the same eclipse (within 30 days)
                if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                    record("P1", f"{ecl_name} {loc_name}", True, "no eclipse visible")
                    continue

                if se_tret[0] == 0.0 or le_tret[0] == 0.0:
                    record(
                        "P1",
                        f"{ecl_name} {loc_name}",
                        False,
                        f"visibility mismatch SE_max={se_tret[0]:.1f} LE_max={le_tret[0]:.1f}",
                    )
                    continue

                if abs(se_tret[0] - le_tret[0]) > 30:
                    record(
                        "P1",
                        f"{ecl_name} {loc_name}",
                        False,
                        f"different eclipses SE={se_tret[0]:.1f} LE={le_tret[0]:.1f}",
                    )
                    continue

                compare_lun_times(
                    "P1",
                    f"{ecl_name} {loc_name}",
                    se_tret,
                    le_tret,
                    max_tol_s=120,
                    phase_tol_s=600,
                )

            except Exception as e:
                errors += 1
                results["errors"].append(f"P1 {ecl_name} {loc_name}: {e}")


def phase2():
    """lun_eclipse_how at known eclipse maximum times."""
    global errors
    print("\n=== P2: lun_eclipse_how at maximum ===")

    test_cases = {
        "Total_2022Nov08_London": (2459890.0, 51.5, -0.1),
        "Total_2025Mar14_NewYork": (2460745.0, 40.7, -74.0),
        "Partial_2024Sep18_Tokyo": (2460572.0, 35.7, 139.7),
        "Total_2019Jan21_SaoPaulo": (2458505.0, -23.5, -46.6),
        "Penumbral_2024Mar25_Sydney": (2460394.0, -33.9, 151.2),
    }

    for name, (jd_before, lat, lon) in test_cases.items():
        try:
            # Get eclipse max time from SE
            geopos = [lon, lat, 0.0]
            se_when = swe.lun_eclipse_when_loc(jd_before, geopos, FLG_SWIEPH, False)
            jd_max = se_when[1][0]

            if jd_max == 0.0:
                record("P2", f"{name}_how", True, "no eclipse visible")
                continue

            result, err = safe_lun_how(jd_max, lat, lon)
            if result is None:
                errors += 1
                results["errors"].append(f"P2 {name}: {err}")
                continue

            se_flag, se_attr, le_flag, le_attr = result
            compare_lun_attrs("P2", f"{name}_how", se_attr, le_attr)

        except Exception as e:
            errors += 1
            results["errors"].append(f"P2 {name}: {e}")


def phase3():
    """Eclipse attributes comparison."""
    global errors
    print("\n=== P3: Eclipse attributes ===")

    test_locs = {
        k: v
        for k, v in LOCATIONS.items()
        if k in ("London", "NewYork", "Tokyo", "CapeTown")
    }

    # Use a few eclipses
    test_eclipses = {
        "Total_2022Nov08": 2459890.0,
        "Partial_2024Sep18": 2460572.0,
        "Total_2025Mar14": 2460745.0,
    }

    for ecl_name, jd_before in test_eclipses.items():
        for loc_name, (lat, lon) in test_locs.items():
            try:
                result, err = safe_lun_when_loc(jd_before, lat, lon)
                if result is None:
                    continue

                se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

                if se_tret[0] == 0.0 or le_tret[0] == 0.0:
                    continue

                if abs(se_tret[0] - le_tret[0]) > 30:
                    continue

                compare_lun_attrs("P3", f"{ecl_name} {loc_name}_attr", se_attr, le_attr)
            except Exception as e:
                errors += 1
                results["errors"].append(f"P3 {ecl_name} {loc_name}: {e}")


def phase4():
    """Moonrise/moonset eclipses."""
    global errors
    print("\n=== P4: Moonrise/moonset eclipses ===")

    # Locations where Moon may rise or set during eclipse
    horizon_cases = {
        "2022Nov08_London": (2459890.0, 51.5, -0.1),
        "2022Nov08_Delhi": (2459890.0, 28.6, 77.2),
        "2025Mar14_Tokyo": (2460745.0, 35.7, 139.7),
        "2024Sep18_LA": (2460572.0, 34.1, -118.2),
        "2019Jan21_Tokyo": (2458505.0, 35.7, 139.7),
        "2021Nov19_London": (2459535.0, 51.5, -0.1),
    }

    for name, (jd_before, lat, lon) in horizon_cases.items():
        try:
            result, err = safe_lun_when_loc(jd_before, lat, lon)
            if result is None:
                errors += 1
                results["errors"].append(f"P4 {name}: {err}")
                continue

            se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

            if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                record("P4", f"{name}", True, "no eclipse visible")
                continue

            if se_tret[0] == 0.0 or le_tret[0] == 0.0:
                record(
                    "P4",
                    f"{name}",
                    False,
                    f"visibility mismatch SE={se_tret[0]:.1f} LE={le_tret[0]:.1f}",
                )
                continue

            if abs(se_tret[0] - le_tret[0]) > 30:
                continue

            # Wider tolerance for horizon events
            compare_lun_times(
                "P4", f"{name}_times", se_tret, le_tret, max_tol_s=180, phase_tol_s=600
            )

        except Exception as e:
            errors += 1
            results["errors"].append(f"P4 {name}: {e}")


def phase5():
    """Southern hemisphere and equatorial eclipses."""
    global errors
    print("\n=== P5: Southern/equatorial locations ===")

    south_locs = {
        "Sydney": (-33.9, 151.2),
        "CapeTown": (-33.9, 18.4),
        "SaoPaulo": (-23.5, -46.6),
        "BuenosAires": (-34.6, -58.4),
        "Auckland": (-36.8, 174.8),
        "Equator0": (0.0, 0.0),
        "Singapore": (1.3, 103.8),
        "Nairobi": (-1.3, 36.8),
    }

    test_eclipses = {
        "Total_2022May16": 2459710.0,
        "Total_2025Sep07": 2460923.0,
        "Total_2018Jan31": 2458150.0,
    }

    for ecl_name, jd_before in test_eclipses.items():
        for loc_name, (lat, lon) in south_locs.items():
            try:
                result, err = safe_lun_when_loc(jd_before, lat, lon)
                if result is None:
                    continue

                se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

                if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                    record("P5", f"{ecl_name} {loc_name}", True, "no eclipse visible")
                    continue

                if se_tret[0] == 0.0 or le_tret[0] == 0.0:
                    record(
                        "P5",
                        f"{ecl_name} {loc_name}",
                        False,
                        f"vis mismatch SE={se_tret[0]:.1f} LE={le_tret[0]:.1f}",
                    )
                    continue

                if abs(se_tret[0] - le_tret[0]) > 30:
                    continue

                compare_lun_times(
                    "P5",
                    f"{ecl_name} {loc_name}",
                    se_tret,
                    le_tret,
                    max_tol_s=120,
                    phase_tol_s=600,
                )

            except Exception as e:
                errors += 1
                results["errors"].append(f"P5 {ecl_name} {loc_name}: {e}")


def phase6():
    """Penumbral-only eclipses."""
    global errors
    print("\n=== P6: Penumbral eclipses ===")

    penumbral = {
        "Penumbral_2024Mar25": 2460394.0,
        "Penumbral_2023May05": 2460069.0,
    }

    test_locs = {
        "London": (51.5, -0.1),
        "Tokyo": (35.7, 139.7),
        "Sydney": (-33.9, 151.2),
        "NewYork": (40.7, -74.0),
        "Delhi": (28.6, 77.2),
    }

    for ecl_name, jd_before in penumbral.items():
        for loc_name, (lat, lon) in test_locs.items():
            try:
                result, err = safe_lun_when_loc(jd_before, lat, lon)
                if result is None:
                    continue

                se_flag, se_tret, se_attr, le_flag, le_tret, le_attr = result

                if se_tret[0] == 0.0 and le_tret[0] == 0.0:
                    record("P6", f"{ecl_name} {loc_name}", True, "no eclipse visible")
                    continue

                if se_tret[0] == 0.0 or le_tret[0] == 0.0:
                    record(
                        "P6",
                        f"{ecl_name} {loc_name}",
                        False,
                        f"vis mismatch SE={se_tret[0]:.1f} LE={le_tret[0]:.1f}",
                    )
                    continue

                if abs(se_tret[0] - le_tret[0]) > 30:
                    continue

                compare_lun_times(
                    "P6",
                    f"{ecl_name} {loc_name}_times",
                    se_tret,
                    le_tret,
                    max_tol_s=120,
                    phase_tol_s=600,
                )
                compare_lun_attrs("P6", f"{ecl_name} {loc_name}_attr", se_attr, le_attr)

            except Exception as e:
                errors += 1
                results["errors"].append(f"P6 {ecl_name} {loc_name}: {e}")


def main():
    print("=" * 70)
    print("ROUND 49: Lunar Eclipse Local Multi-Location Sweep")
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
    print(f"ROUND 49 FINAL: {passed}/{total} passed ({pct:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print("=" * 70)

    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:50]:
            print(f)

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:10]:
            print(e)

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
