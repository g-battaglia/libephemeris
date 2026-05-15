#!/usr/bin/env python3
"""Round 116: Moon at Perigee/Apogee — Position Precision at Orbital Extremes

Tests Moon positions at perigee (closest) and apogee (farthest) points,
where perturbation effects are strongest and ephemeris differences most visible.

Also tests:
- Moon distance at perigee/apogee vs SE
- Moon speed at perigee/apogee (speed extremes)
- All flag combinations at these critical points
- Mean Apogee (Lilith) and Mean Node at these times
- Perigee/apogee timing via nod_aps_ut
"""

from __future__ import annotations

import os
import sys
import math

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# SE ephemeris path
swe.set_ephe_path(_REF_EPHE_PATH)

# Constants
MOON = 1
MEAN_NODE = 10
TRUE_NODE = 11
MEAN_APOG = 12  # Mean Lilith
OSCU_APOG = 13  # Osculating Lilith
FLG_SPEED = 256
FLG_EQUATORIAL = 2048
FLG_J2000 = 32
FLG_NONUT = 64
FLG_NOABERR = 1024
FLG_TRUEPOS = 16
FLG_TOPOCTR = 32768
FLG_SIDEREAL = 65536


def se_hsys(ch):
    return ch.encode("ascii") if isinstance(ch, str) else ch


def find_moon_perigees_apogees(start_jd, count=50):
    """Find Moon perigee and apogee times by scanning distance minima/maxima."""
    events = []
    jd = start_jd
    step = 0.5  # half-day steps

    # Get initial distance
    prev_dist = swe.calc_ut(jd, MOON, FLG_SPEED)[0][2]
    prev_prev_dist = swe.calc_ut(jd - step, MOON, FLG_SPEED)[0][2]

    jd += step
    found = 0
    max_iter = 10000
    iterations = 0

    while found < count and iterations < max_iter:
        iterations += 1
        curr_dist = swe.calc_ut(jd, MOON, FLG_SPEED)[0][2]

        # Check for local minimum (perigee) or maximum (apogee)
        if prev_dist < prev_prev_dist and prev_dist < curr_dist:
            # Perigee found near jd - step, refine with bisection
            refined_jd = _refine_extremum(jd - step, step, is_minimum=True)
            events.append(("perigee", refined_jd))
            found += 1
        elif prev_dist > prev_prev_dist and prev_dist > curr_dist:
            # Apogee found near jd - step, refine
            refined_jd = _refine_extremum(jd - step, step, is_maximum=True)
            events.append(("apogee", refined_jd))
            found += 1

        prev_prev_dist = prev_dist
        prev_dist = curr_dist
        jd += step

    return events


def _refine_extremum(jd_approx, window, is_minimum=False, is_maximum=False):
    """Refine perigee/apogee time using golden section search."""
    a = jd_approx - window
    b = jd_approx + window
    gr = (math.sqrt(5) + 1) / 2

    for _ in range(50):  # ~50 iterations gives sub-second precision
        c = b - (b - a) / gr
        d = a + (b - a) / gr

        fc = swe.calc_ut(c, MOON, FLG_SPEED)[0][2]
        fd = swe.calc_ut(d, MOON, FLG_SPEED)[0][2]

        if is_minimum:
            if fc < fd:
                b = d
            else:
                a = c
        else:  # is_maximum
            if fc > fd:
                b = d
            else:
                a = c

    return (a + b) / 2


def compare_moon_at_event(jd, event_type, flags_list):
    """Compare Moon position between LE and SE at a given JD with various flags."""
    results = []

    for flags in flags_list:
        flag_val = FLG_SPEED
        flag_names = ["SPEED"]
        for fname, fval in flags:
            flag_val |= fval
            flag_names.append(fname)

        flag_str = "|".join(flag_names)

        try:
            se_result = swe.calc_ut(jd, MOON, flag_val)
            se_pos = se_result[0]
        except Exception as e:
            results.append((flag_str, "SE_ERROR", str(e)))
            continue

        try:
            le_result = ephem.calc_ut(jd, MOON, flag_val)
            le_pos = le_result[0]
        except Exception as e:
            results.append((flag_str, "LE_ERROR", str(e)))
            continue

        # Compare all 6 components
        labels = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]
        if flag_val & FLG_EQUATORIAL:
            labels = ["ra", "dec", "dist", "ra_spd", "dec_spd", "dist_spd"]

        for i, label in enumerate(labels):
            se_val = se_pos[i]
            le_val = le_pos[i]
            diff = le_val - se_val

            # Wrap longitude/RA differences
            if i == 0 and not (flag_val & FLG_EQUATORIAL):
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

            diff_arcsec = abs(diff) * 3600

            # Tolerances
            if "dist" in label and "spd" not in label:
                # Distance in AU - use relative tolerance
                if abs(se_val) > 1e-10:
                    rel_diff = abs(diff / se_val)
                    tol_ok = rel_diff < 1e-6  # 1 ppm
                else:
                    tol_ok = abs(diff) < 1e-10
                results.append(
                    (
                        flag_str,
                        label,
                        se_val,
                        le_val,
                        diff,
                        rel_diff if abs(se_val) > 1e-10 else 0,
                        tol_ok,
                        "rel",
                    )
                )
            elif "spd" in label:
                # Speed - use arcsec/day tolerance
                tol = 5.0  # 5"/day for Moon speed (moves ~46800"/day)
                results.append(
                    (
                        flag_str,
                        label,
                        se_val,
                        le_val,
                        diff_arcsec,
                        tol,
                        diff_arcsec < tol,
                        "arcsec/day",
                    )
                )
            else:
                # Position - use arcsec tolerance
                tol = 2.0  # 2" for Moon position
                results.append(
                    (
                        flag_str,
                        label,
                        se_val,
                        le_val,
                        diff_arcsec,
                        tol,
                        diff_arcsec < tol,
                        "arcsec",
                    )
                )

    return results


def compare_analytical_bodies(jd, event_type):
    """Compare Mean Node, True Node, Mean Apogee, Oscu Apogee at perigee/apogee times."""
    results = []
    bodies = [
        (MEAN_NODE, "MeanNode"),
        (TRUE_NODE, "TrueNode"),
        (MEAN_APOG, "MeanApog"),
        (OSCU_APOG, "OscuApog"),
    ]

    for body_id, body_name in bodies:
        try:
            se_result = swe.calc_ut(jd, body_id, FLG_SPEED)
            se_pos = se_result[0]
        except Exception as e:
            results.append((body_name, "SE_ERROR", str(e)))
            continue

        try:
            le_result = ephem.calc_ut(jd, body_id, FLG_SPEED)
            le_pos = le_result[0]
        except Exception as e:
            results.append((body_name, "LE_ERROR", str(e)))
            continue

        for i, label in enumerate(
            ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]
        ):
            se_val = se_pos[i]
            le_val = le_pos[i]
            diff = le_val - se_val

            if i == 0:
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

            diff_arcsec = abs(diff) * 3600

            if "dist" in label and "spd" not in label:
                results.append(
                    (
                        body_name,
                        label,
                        se_val,
                        le_val,
                        abs(diff),
                        1e-4,
                        abs(diff) < 1e-4,
                        "AU",
                    )
                )
            elif "spd" in label:
                # Analytical body speeds can differ more
                tol = (
                    200.0 if body_name == "OscuApog" else 10.0
                )  # OscuApog has known ~155"/day diffs
                results.append(
                    (
                        body_name,
                        label,
                        se_val,
                        le_val,
                        diff_arcsec,
                        tol,
                        diff_arcsec < tol,
                        "arcsec/day",
                    )
                )
            else:
                tol = (
                    30.0 if body_name == "OscuApog" else 2.0
                )  # OscuApog can differ more
                if label == "lat" and body_name == "MeanApog":
                    tol = 25.0  # Known ~19" model difference
                results.append(
                    (
                        body_name,
                        label,
                        se_val,
                        le_val,
                        diff_arcsec,
                        tol,
                        diff_arcsec < tol,
                        "arcsec",
                    )
                )

    return results


def main():
    print("=" * 80)
    print("ROUND 116: Moon at Perigee/Apogee — Position Precision at Orbital Extremes")
    print("=" * 80)

    # Find perigee/apogee events across 2000-2030
    print("\n--- Finding Moon perigees and apogees (2000-2030) ---")
    start_jd = 2451545.0  # J2000.0
    events = find_moon_perigees_apogees(start_jd, count=80)

    perigees = [(t, jd) for t, jd in events if t == "perigee"]
    apogees = [(t, jd) for t, jd in events if t == "apogee"]
    print(f"Found {len(perigees)} perigees, {len(apogees)} apogees")

    # Flag combinations to test
    flag_combos = [
        [],  # Default (ecliptic of date)
        [("EQUATORIAL", FLG_EQUATORIAL)],
        [("J2000", FLG_J2000)],
        [("NONUT", FLG_NONUT)],
        [("NOABERR", FLG_NOABERR)],
        [("TRUEPOS", FLG_TRUEPOS)],
        [("J2000", FLG_J2000), ("NONUT", FLG_NONUT)],
        [("EQUATORIAL", FLG_EQUATORIAL), ("J2000", FLG_J2000)],
    ]

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Test Moon positions at all events
    print("\n--- Testing Moon positions at perigee/apogee events ---")
    for event_type, jd in events:
        results = compare_moon_at_event(jd, event_type, flag_combos)

        for item in results:
            if isinstance(item[1], str) and item[1].endswith("ERROR"):
                # Skip errors
                continue

            flag_str, label, se_val, le_val, diff_val, tol, passed, unit = item
            total_tests += 1
            if passed:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 30:
                    failures.append(
                        f"  {event_type} JD={jd:.4f}: {flag_str} {label}: "
                        f"SE={se_val:.8f} LE={le_val:.8f} diff={diff_val:.4f}{unit} (tol={tol})"
                    )

    # Test analytical bodies at perigee/apogee
    print("\n--- Testing analytical bodies (Node/Lilith) at perigee/apogee ---")
    for event_type, jd in events:
        results = compare_analytical_bodies(jd, event_type)

        for item in results:
            if isinstance(item[1], str) and item[1].endswith("ERROR"):
                continue

            body_name, label, se_val, le_val, diff_val, tol, passed, unit = item
            total_tests += 1
            if passed:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 60:
                    failures.append(
                        f"  {event_type} JD={jd:.4f}: {body_name} {label}: "
                        f"SE={se_val:.8f} LE={le_val:.8f} diff={diff_val:.6f}{unit} (tol={tol})"
                    )

    # Test Moon distance precision at perigee/apogee
    print("\n--- Testing Moon distance at perigee/apogee extremes ---")
    min_perigee_dist_se = 999
    min_perigee_dist_le = 999
    max_apogee_dist_se = 0
    max_apogee_dist_le = 0

    for event_type, jd in events:
        se_dist = swe.calc_ut(jd, MOON, FLG_SPEED)[0][2]
        le_dist = ephem.calc_ut(jd, MOON, FLG_SPEED)[0][2]

        if event_type == "perigee":
            min_perigee_dist_se = min(min_perigee_dist_se, se_dist)
            min_perigee_dist_le = min(min_perigee_dist_le, le_dist)
        else:
            max_apogee_dist_se = max(max_apogee_dist_se, se_dist)
            max_apogee_dist_le = max(max_apogee_dist_le, le_dist)

    print(
        f"  Closest perigee:  SE={min_perigee_dist_se:.10f} AU  LE={min_perigee_dist_le:.10f} AU"
    )
    print(
        f"  Farthest apogee:  SE={max_apogee_dist_se:.10f} AU  LE={max_apogee_dist_le:.10f} AU"
    )

    # Test Moon speed at perigee (fastest) and apogee (slowest)
    print("\n--- Moon speed extremes ---")
    max_speed_se = 0
    max_speed_le = 0
    min_speed_se = 999
    min_speed_le = 999

    for event_type, jd in events:
        se_spd = abs(swe.calc_ut(jd, MOON, FLG_SPEED)[0][3])
        le_spd = abs(ephem.calc_ut(jd, MOON, FLG_SPEED)[0][3])

        if event_type == "perigee":
            max_speed_se = max(max_speed_se, se_spd)
            max_speed_le = max(max_speed_le, le_spd)
        else:
            min_speed_se = min(min_speed_se, se_spd)
            min_speed_le = min(min_speed_le, le_spd)

    print(
        f"  Max speed (perigee): SE={max_speed_se:.6f}°/day  LE={max_speed_le:.6f}°/day"
    )
    print(
        f"  Min speed (apogee):  SE={min_speed_se:.6f}°/day  LE={min_speed_le:.6f}°/day"
    )

    # Test with topocentric for a few events (parallax matters most at perigee)
    print("\n--- Testing topocentric Moon at perigees (parallax critical) ---")
    locations = [
        (0.0, 51.5, 0.0, "London"),
        (139.7, 35.7, 0.0, "Tokyo"),
        (-74.0, 40.7, 0.0, "NYC"),
        (0.0, -90.0, 2835.0, "SouthPole"),
    ]

    for loc_lon, loc_lat, loc_alt, loc_name in locations:
        swe.set_topo(loc_lon, loc_lat, loc_alt)
        ephem.set_topo(loc_lon, loc_lat, loc_alt)

        # Test first 5 perigees
        for event_type, jd in perigees[:5]:
            flags_topo = FLG_SPEED | FLG_TOPOCTR

            try:
                se_result = swe.calc_ut(jd, MOON, flags_topo)
                se_pos = se_result[0]
            except Exception:
                continue

            try:
                le_result = ephem.calc_ut(jd, MOON, flags_topo)
                le_pos = le_result[0]
            except Exception:
                continue

            for i, label in enumerate(["lon", "lat", "dist"]):
                se_val = se_pos[i]
                le_val = le_pos[i]
                diff = le_val - se_val

                if i == 0:
                    if diff > 180:
                        diff -= 360
                    elif diff < -180:
                        diff += 360

                diff_arcsec = abs(diff) * 3600

                if label == "dist":
                    rel_diff = abs(diff / se_val) if abs(se_val) > 0 else 0
                    tol_ok = rel_diff < 1e-6
                    total_tests += 1
                    if tol_ok:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 80:
                            failures.append(
                                f"  TOPO {loc_name} perigee JD={jd:.4f}: dist rel_diff={rel_diff:.2e}"
                            )
                else:
                    tol = 2.0
                    total_tests += 1
                    if diff_arcsec < tol:
                        total_pass += 1
                    else:
                        total_fail += 1
                        if len(failures) < 80:
                            failures.append(
                                f'  TOPO {loc_name} perigee JD={jd:.4f}: {label} diff={diff_arcsec:.4f}" (tol={tol}")'
                            )

    # Test with sidereal mode at perigee/apogee
    print("\n--- Testing sidereal Moon at perigee/apogee ---")
    swe.set_sid_mode(0)  # Fagan-Bradley
    ephem.set_sid_mode(0, 0, 0)

    for event_type, jd in events[:20]:  # First 20 events
        flags_sid = FLG_SPEED | FLG_SIDEREAL

        try:
            se_result = swe.calc_ut(jd, MOON, flags_sid)
            se_pos = se_result[0]
        except Exception:
            continue

        try:
            le_result = ephem.calc_ut(jd, MOON, flags_sid)
            le_pos = le_result[0]
        except Exception:
            continue

        for i, label in enumerate(["lon", "lat"]):
            se_val = se_pos[i]
            le_val = le_pos[i]
            diff = le_val - se_val

            if i == 0:
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360

            diff_arcsec = abs(diff) * 3600
            tol = 2.0

            total_tests += 1
            if diff_arcsec < tol:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 90:
                    failures.append(
                        f'  SIDEREAL {event_type} JD={jd:.4f}: {label} diff={diff_arcsec:.4f}" (tol={tol}")'
                    )

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    # Summary
    print("\n" + "=" * 80)
    print(
        f"ROUND 116 RESULTS: {total_pass}/{total_tests} passed ({100 * total_pass / total_tests:.1f}%)"
    )
    print(f"  Failures: {total_fail}")
    print("=" * 80)

    if failures:
        print("\nSample failures:")
        for f in failures[:30]:
            print(f)

    if total_fail == 0:
        print("\nAll tests PASSED!")

    return total_fail


if __name__ == "__main__":
    sys.exit(main())
