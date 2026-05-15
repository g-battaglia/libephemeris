#!/usr/bin/env python3
"""Round 152: Eclipse magnitude at contact times.

Compare sol_eclipse_how at known eclipse maximum times to verify
eclipse attributes (magnitude, obscuration, etc.) match between
libephemeris and pyswisseph.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256

# Known eclipse dates and locations
# (jd_approx, lat, lon, alt, description)
eclipse_cases = []

# Find eclipses using sol_eclipse_when_glob first
jd = 2451545.0  # J2000
for i in range(30):  # Find 30 eclipses
    try:
        se_r = swe.sol_eclipse_when_glob(jd, FLG_SPEED, 0, False)
        t_max = se_r[1][0]
        if t_max > 0:
            # Test at several locations
            for lat, lon, loc_name in [
                (0.0, 0.0, "equator"),
                (45.0, 0.0, "45N"),
                (-30.0, 0.0, "30S"),
                (52.0, 13.4, "Berlin"),
                (35.7, 139.7, "Tokyo"),
            ]:
                eclipse_cases.append(
                    (t_max, lat, lon, 0.0, f"Eclipse#{i + 1} {loc_name}")
                )
            jd = t_max + 20  # next eclipse
        else:
            break
    except:
        jd += 180
        continue

# Attribute labels
ATTR_LABELS = [
    "ecl_magnitude",
    "surface_ratio",
    "disc_fraction",
    "diam_ratio",
    "moon_sun_dist",
    "saros_nr",
    "saros_member",
    "attr7",
    "attr8",
    "attr9",
]

# Tolerances
TOL_MAG = 0.01  # magnitude
TOL_RATIO = 0.01  # surface/disc ratios
TOL_DIAM = 0.001  # diameter ratio
TOL_DIST = 0.01  # degrees for moon-sun distance

passed = failed = errors = total = 0
failures = []

print(f"Round 152: Eclipse Magnitude at Contact Times")
print(f"Testing {len(eclipse_cases)} cases")
print("=" * 90)

for t_max, lat, lon, alt, desc in eclipse_cases:
    try:
        se_r = swe.sol_eclipse_how(t_max, [lon, lat, alt], FLG_SPEED)
        le_r = ephem.sol_eclipse_how(t_max, FLG_SPEED, [lon, lat, alt])

        se_retflag = se_r[0]
        le_retflag = le_r[0]
        se_attr = se_r[1]
        le_attr = le_r[1]

        # Compare first 5 attributes (the meaningful ones)
        for i in range(min(5, len(se_attr), len(le_attr))):
            total += 1
            se_val = se_attr[i]
            le_val = le_attr[i]

            # Skip if both zero (no eclipse visible)
            if se_val == 0.0 and le_val == 0.0:
                passed += 1
                continue

            if i == 0:  # magnitude
                diff = abs(le_val - se_val)
                tol = TOL_MAG
            elif i <= 2:  # surface/disc ratio
                diff = abs(le_val - se_val)
                tol = TOL_RATIO
            elif i == 3:  # diam ratio
                diff = abs(le_val - se_val)
                tol = TOL_DIAM
            else:  # moon-sun dist
                diff = abs(le_val - se_val)
                tol = TOL_DIST

            if diff <= tol:
                passed += 1
            else:
                failed += 1
                aname = ATTR_LABELS[i] if i < len(ATTR_LABELS) else f"attr{i}"
                msg = f"  FAIL {desc} {aname}: SE={se_val:.8f} LE={le_val:.8f} diff={diff:.6f}"
                failures.append(msg)
                if len(failures) <= 20:
                    print(msg)

    except Exception as e:
        errors += 1
        if errors <= 3:
            print(f"  ERROR {desc}: {e}")

print()
print("=" * 90)
if total > 0:
    print(
        f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
    )
else:
    print(f"No tests run. {errors} errors.")
if failures:
    print(f"\n{len(failures)} failures total")
else:
    print("\nAll tests passed!")
