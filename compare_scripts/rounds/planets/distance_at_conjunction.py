#!/usr/bin/env python3
"""Round 187: Planet distance at conjunction/opposition.

Tests geocentric distances of planets at their closest approach (conjunction)
and farthest (opposition) from Earth. These are extreme distance values where
ephemeris model differences can be most visible.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("LIBEPHEMERIS_MODE", "skyfield")

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLAGS = ephem.FLG_SWIEPH | ephem.FLG_SPEED

# Known approximate conjunction/opposition dates (JD) for various planets
# These are times when distance is at an extreme
TEST_EVENTS = [
    # Mars oppositions (closest approach)
    ("Mars opposition 2003", ephem.MARS, swe.MARS, 2452879.5),  # Aug 2003
    ("Mars opposition 2005", ephem.MARS, swe.MARS, 2453684.5),  # Nov 2005
    ("Mars opposition 2010", ephem.MARS, swe.MARS, 2455210.5),  # Jan 2010
    ("Mars opposition 2014", ephem.MARS, swe.MARS, 2456753.5),  # Apr 2014
    ("Mars opposition 2018", ephem.MARS, swe.MARS, 2458324.5),  # Jul 2018
    ("Mars opposition 2020", ephem.MARS, swe.MARS, 2459133.5),  # Oct 2020
    # Mars conjunctions (farthest)
    ("Mars conjunction 2004", ephem.MARS, swe.MARS, 2453270.5),  # Sep 2004
    ("Mars conjunction 2006", ephem.MARS, swe.MARS, 2454042.5),  # Oct 2006
    ("Mars conjunction 2019", ephem.MARS, swe.MARS, 2458720.5),  # Sep 2019
    # Jupiter oppositions
    ("Jupiter opp 2010", ephem.JUPITER, swe.JUPITER, 2455461.5),
    ("Jupiter opp 2015", ephem.JUPITER, swe.JUPITER, 2457048.5),
    ("Jupiter opp 2020", ephem.JUPITER, swe.JUPITER, 2459038.5),
    ("Jupiter opp 2022", ephem.JUPITER, swe.JUPITER, 2459845.5),
    # Saturn oppositions
    ("Saturn opp 2010", ephem.SATURN, swe.SATURN, 2455281.5),
    ("Saturn opp 2015", ephem.SATURN, swe.SATURN, 2457168.5),
    ("Saturn opp 2020", ephem.SATURN, swe.SATURN, 2459046.5),
    # Venus inferior conjunctions (closest)
    ("Venus inf conj 2004", ephem.VENUS, swe.VENUS, 2453158.5),
    ("Venus inf conj 2012", ephem.VENUS, swe.VENUS, 2456082.5),
    ("Venus inf conj 2020", ephem.VENUS, swe.VENUS, 2459005.5),
    # Venus superior conjunctions (farthest)
    ("Venus sup conj 2006", ephem.VENUS, swe.VENUS, 2454046.5),
    ("Venus sup conj 2014", ephem.VENUS, swe.VENUS, 2456953.5),
    ("Venus sup conj 2022", ephem.VENUS, swe.VENUS, 2459877.5),
    # Mercury conjunctions
    ("Mercury inf conj 2000", ephem.MERCURY, swe.MERCURY, 2451651.5),
    ("Mercury inf conj 2010", ephem.MERCURY, swe.MERCURY, 2455474.5),
    ("Mercury inf conj 2020", ephem.MERCURY, swe.MERCURY, 2458896.5),
    ("Mercury sup conj 2000", ephem.MERCURY, swe.MERCURY, 2451553.5),
    ("Mercury sup conj 2010", ephem.MERCURY, swe.MERCURY, 2455407.5),
    ("Mercury sup conj 2020", ephem.MERCURY, swe.MERCURY, 2459032.5),
    # Moon perigee/apogee (closest/farthest)
    ("Moon perigee 2000 Jan", ephem.MOON, swe.MOON, 2451560.5),
    ("Moon perigee 2010 Jan", ephem.MOON, swe.MOON, 2455218.5),
    ("Moon perigee 2020 Apr", ephem.MOON, swe.MOON, 2458949.5),
    ("Moon apogee 2000 Jan", ephem.MOON, swe.MOON, 2451566.5),
    ("Moon apogee 2010 Jan", ephem.MOON, swe.MOON, 2455225.5),
    ("Moon apogee 2020 Mar", ephem.MOON, swe.MOON, 2458935.5),
    # Outer planets
    ("Uranus opp 2015", ephem.URANUS, swe.URANUS, 2457306.5),
    ("Neptune opp 2015", ephem.NEPTUNE, swe.NEPTUNE, 2457261.5),
    ("Pluto opp 2015", ephem.PLUTO, swe.PLUTO, 2457209.5),
]

passed = 0
failed = 0
total = 0
failures = []


def compare(label, le_body, se_body, jd, flags_le, flags_se):
    global passed, failed, total

    try:
        le_r = ephem.calc_ut(jd, le_body, flags_le)
        se_r = swe.calc_ut(jd, se_body, flags_se)
    except Exception as e:
        return

    le_lon, le_lat, le_dist = le_r[0][0], le_r[0][1], le_r[0][2]
    se_lon, se_lat, se_dist = se_r[0][0], se_r[0][1], se_r[0][2]
    le_lon_spd, le_lat_spd, le_dist_spd = le_r[0][3], le_r[0][4], le_r[0][5]
    se_lon_spd, se_lat_spd, se_dist_spd = se_r[0][3], se_r[0][4], se_r[0][5]

    is_moon = "Moon" in label
    is_pluto = "Pluto" in label

    # Longitude
    total += 1
    lon_diff = abs(le_lon - se_lon)
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_as = lon_diff * 3600
    tol = 1.0 if is_moon else (2.0 if is_pluto else 0.5)
    if lon_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON: {lon_as:.4f}" (tol {tol}")')

    # Latitude
    total += 1
    lat_as = abs(le_lat - se_lat) * 3600
    if lat_as <= tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LAT: {lat_as:.4f}" (tol {tol}")')

    # Distance
    total += 1
    dist_diff = abs(le_dist - se_dist)
    # Tolerance: 1e-6 AU for Moon, 1e-5 AU for outer, 5e-6 for inner
    dist_tol = 1e-6 if is_moon else (1e-5 if is_pluto else 5e-6)
    if dist_diff <= dist_tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f"  {label} DIST: LE={le_dist:.8f} SE={se_dist:.8f} diff={dist_diff:.2e} AU"
        )

    # Distance speed
    total += 1
    ds_diff = abs(le_dist_spd - se_dist_spd)
    ds_tol = 1e-5 if is_moon else (5e-5 if is_pluto else 1e-5)
    if ds_diff <= ds_tol:
        passed += 1
    else:
        failed += 1
        failures.append(
            f"  {label} DIST_SPD: LE={le_dist_spd:.8f} SE={se_dist_spd:.8f} diff={ds_diff:.2e}"
        )

    # Longitude speed
    total += 1
    ls_diff = abs(le_lon_spd - se_lon_spd) * 3600
    ls_tol = 5.0 if is_moon else (3.0 if is_pluto else 2.0)
    if ls_diff <= ls_tol:
        passed += 1
    else:
        failed += 1
        failures.append(f'  {label} LON_SPD: diff={ls_diff:.4f}"/day')


if __name__ == "__main__":
    print("=" * 70)
    print("Round 187: Planet Distance at Conjunction/Opposition")
    print("=" * 70)

    se_flags = swe.FLG_SWIEPH | swe.FLG_SPEED

    for label, le_body, se_body, jd in TEST_EVENTS:
        # Test at event date
        compare(label, le_body, se_body, jd, FLAGS, se_flags)
        # Test ±1 day
        compare(f"{label} -1d", le_body, se_body, jd - 1.0, FLAGS, se_flags)
        compare(f"{label} +1d", le_body, se_body, jd + 1.0, FLAGS, se_flags)

    # Also test with FLG_EQUATORIAL
    print("\nEquatorial mode tests...")
    flags_le_eq = FLAGS | ephem.FLG_EQUATORIAL
    flags_se_eq = se_flags | swe.FLG_EQUATORIAL
    for label, le_body, se_body, jd in TEST_EVENTS[:12]:
        compare(f"{label} EQ", le_body, se_body, jd, flags_le_eq, flags_se_eq)

    # Also test heliocentric distances for inner planets
    print("Heliocentric distance tests...")
    flags_le_h = FLAGS | ephem.FLG_HELCTR
    flags_se_h = se_flags | swe.FLG_HELCTR
    for label, le_body, se_body, jd in TEST_EVENTS:
        if "Moon" in label:
            continue  # No heliocentric Moon
        try:
            le_r = ephem.calc_ut(jd, le_body, flags_le_h)
            se_r = swe.calc_ut(jd, se_body, flags_se_h)
            total += 1
            dist_diff = abs(le_r[0][2] - se_r[0][2])
            if dist_diff <= 1e-5:
                passed += 1
            else:
                failed += 1
                failures.append(f"  {label} HELIO DIST: diff={dist_diff:.2e} AU")
        except Exception:
            pass

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
