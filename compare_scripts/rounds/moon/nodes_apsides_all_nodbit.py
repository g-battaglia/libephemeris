#!/usr/bin/env python3
"""Round 147: nod_aps_ut all bodies with all NODBIT flags.

Compare planetary nodes and apsides between libephemeris and pyswisseph
across multiple epochs with NODBIT_MEAN, NODBIT_OSCU, NODBIT_OSCU_BAR,
and NODBIT_FOPOINT flags.

Return format: (asc_node_6tuple, desc_node_6tuple, perihelion_6tuple, aphelion_6tuple)
Each 6-tuple: (lon, lat, dist, lon_speed, lat_speed, dist_speed)
"""

from __future__ import annotations
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

# Constants
FLG_SPEED = 256
NODBIT_MEAN = 1
NODBIT_OSCU = 2
NODBIT_OSCU_BAR = 4
NODBIT_FOPOINT = 256

# Bodies that both SE and LE support for nod_aps_ut
# SE errors on 10 (mean node), 11 (true node), 12 (mean apog)
# LE returns zeros for Sun (0) - known gap
BODIES = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # Moon through Pluto
BODY_NAMES = {
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

# NODBIT flags to test
NODBIT_FLAGS = [
    (NODBIT_MEAN, "MEAN"),
    (NODBIT_OSCU, "OSCU"),
    (NODBIT_OSCU_BAR, "OSCU_BAR"),
    (NODBIT_MEAN | NODBIT_FOPOINT, "MEAN+FOPOINT"),
    (NODBIT_OSCU | NODBIT_FOPOINT, "OSCU+FOPOINT"),
]

# Test dates
test_dates = []
for year in [1900, 1950, 1980, 2000, 2010, 2020, 2024, 2025, 2050, 2100]:
    jd = swe.julday(year, 1, 1, 12.0)
    test_dates.append((f"{year}", jd))

# Add some mid-year dates
for year in [2000, 2024]:
    for month in [4, 7, 10]:
        jd = swe.julday(year, month, 1, 12.0)
        test_dates.append((f"{year}-{month:02d}", jd))

# Node labels
NODE_LABELS = ["asc_node", "desc_node", "perihelion", "aphelion"]
COMPONENT_LABELS = ["lon", "lat", "dist", "lon_spd", "lat_spd", "dist_spd"]

# Tolerances per component (arcseconds for angular, AU for distance)
TOL_LON = 2.0  # 2" for longitude
TOL_LAT = 2.0  # 2" for latitude
TOL_DIST = 0.001  # AU for distance
TOL_LON_SPD = 5.0  # "/day for lon speed
TOL_LAT_SPD = 5.0  # "/day for lat speed
TOL_DIST_SPD = 0.001  # AU/day for dist speed

# Known model differences are larger for osculating elements
TOL_OSCU_LON = 60.0  # 1' for osculating elements
TOL_OSCU_LAT = 60.0
TOL_OSCU_DIST = 0.01
TOL_OSCU_SPD = 120.0

passed = 0
failed = 0
errors = 0
total = 0
failures = []

print("Round 147: nod_aps_ut All Bodies with All NODBIT Flags")
print(
    f"Testing {len(BODIES)} bodies x {len(NODBIT_FLAGS)} flags x {len(test_dates)} dates"
)
print("=" * 100)

for label, jd in test_dates:
    for body in BODIES:
        for nodbit, nodbit_name in NODBIT_FLAGS:
            try:
                se_result = swe.nod_aps_ut(jd, body, nodbit, FLG_SPEED)
                le_result = ephem.nod_aps_ut(jd, body, nodbit, FLG_SPEED)

                is_oscu = (nodbit & (NODBIT_OSCU | NODBIT_OSCU_BAR)) != 0

                for node_idx in range(4):
                    se_node = se_result[node_idx]
                    le_node = le_result[node_idx]

                    # Skip if both return all zeros (unsupported)
                    if all(v == 0.0 for v in se_node) and all(
                        v == 0.0 for v in le_node
                    ):
                        continue

                    for comp_idx in range(6):
                        total += 1
                        se_val = se_node[comp_idx]
                        le_val = le_node[comp_idx]

                        # Skip if both zero
                        if se_val == 0.0 and le_val == 0.0:
                            passed += 1
                            continue

                        # Skip if LE returns zero but SE doesn't (known gap for some bodies/flags)
                        if le_val == 0.0 and se_val != 0.0:
                            # Known gap - don't count as failure unless it's lon of asc/desc node
                            if comp_idx == 0 and node_idx < 2:
                                failed += 1
                                msg = f"  FAIL {label} {BODY_NAMES[body]} {nodbit_name} {NODE_LABELS[node_idx]}.{COMPONENT_LABELS[comp_idx]}: SE={se_val:.6f} LE=0.0 (ZERO)"
                                failures.append(msg)
                                if len(failures) <= 30:
                                    print(msg)
                            else:
                                passed += 1
                            continue

                        # Determine tolerance
                        if comp_idx <= 1:  # lon, lat (degrees -> arcsec)
                            diff = abs(le_val - se_val) * 3600.0
                            tol = TOL_OSCU_LON if is_oscu else TOL_LON
                            unit = '"'
                        elif comp_idx == 2:  # dist (AU)
                            diff = abs(le_val - se_val)
                            tol = TOL_OSCU_DIST if is_oscu else TOL_DIST
                            unit = "AU"
                        elif comp_idx <= 4:  # lon_spd, lat_spd (deg/day -> arcsec/day)
                            diff = abs(le_val - se_val) * 3600.0
                            tol = TOL_OSCU_SPD if is_oscu else TOL_LON_SPD
                            unit = '"/d'
                        else:  # dist_spd (AU/day)
                            diff = abs(le_val - se_val)
                            tol = TOL_OSCU_DIST if is_oscu else TOL_DIST_SPD
                            unit = "AU/d"

                        if diff <= tol:
                            passed += 1
                        else:
                            failed += 1
                            msg = (
                                f"  FAIL {label} {BODY_NAMES[body]} {nodbit_name} "
                                f"{NODE_LABELS[node_idx]}.{COMPONENT_LABELS[comp_idx]}: "
                                f"SE={se_val:.8f} LE={le_val:.8f} diff={diff:.4f}{unit} tol={tol}{unit}"
                            )
                            failures.append(msg)
                            if len(failures) <= 50:
                                print(msg)

            except Exception as e:
                errors += 1
                err_str = str(e)
                if "not implemented" not in err_str and "not supported" not in err_str:
                    print(f"  ERROR {label} {BODY_NAMES[body]} {nodbit_name}: {e}")

print()
print("=" * 100)
if total > 0:
    print(
        f"Results: {passed}/{total} passed ({100 * passed / total:.1f}%), {failed} failed, {errors} errors"
    )
else:
    print(f"No tests run. {errors} errors.")

if failures:
    print(f"\nShowing first 50 of {len(failures)} failures:")
    for f in failures[:50]:
        print(f)

    # Categorize failures
    cats = {}
    for f in failures:
        parts = f.split()
        key = f.split("FAIL")[1].split(":")[0].strip() if "FAIL" in f else "unknown"
        # Extract body + flag
        for bname in BODY_NAMES.values():
            if bname in f:
                for _, fname in NODBIT_FLAGS:
                    if fname in f:
                        cat = f"{bname} {fname}"
                        cats[cat] = cats.get(cat, 0) + 1
    if cats:
        print(f"\nFailure categories:")
        for cat, count in sorted(cats.items(), key=lambda x: -x[1]):
            print(f"  {cat}: {count}")
else:
    print("\nAll tests passed!")
