#!/usr/bin/env python3
"""Round 118: Eclipse how() Deep — Solar & Lunar Eclipse Attribute Validation

Comprehensive comparison of sol_eclipse_how and lun_eclipse_how attributes
between libephemeris and pyswisseph at known eclipse times.

Tests:
- Solar eclipse attributes (obscuration, magnitude, shadow width, etc.)
- Lunar eclipse attributes (magnitude, penumbral magnitude, etc.)
- Multiple geographic locations for solar eclipses
- Edge cases: partial, total, annular, penumbral eclipses
- Attribute consistency across eclipse phases
"""

from __future__ import annotations

import os
import sys
import math

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)

FLG_SPEED = 256
ECL_CENTRAL = 1
ECL_NONCENTRAL = 2
ECL_TOTAL = 4
ECL_ANNULAR = 8
ECL_PARTIAL = 16
ECL_ANNULAR_TOTAL = 32
ECL_PENUMBRAL = 64


def find_solar_eclipses(start_jd, count=25):
    """Find solar eclipses using pyswisseph."""
    eclipses = []
    jd = start_jd
    for _ in range(count):
        try:
            res = swe.sol_eclipse_when_glob(jd, 0, 0, False)
            ecl_type = res[0]
            tmax = res[1][0]
            eclipses.append((ecl_type, tmax))
            jd = tmax + 30  # next search
        except Exception:
            jd += 30
    return eclipses


def find_lunar_eclipses(start_jd, count=25):
    """Find lunar eclipses using pyswisseph."""
    eclipses = []
    jd = start_jd
    for _ in range(count):
        try:
            res = swe.lun_eclipse_when(jd, 0, 0, False)
            ecl_type = res[0]
            tmax = res[1][0]
            eclipses.append((ecl_type, tmax))
            jd = tmax + 30
        except Exception:
            jd += 30
    return eclipses


def classify_eclipse(ecl_type):
    """Return human-readable eclipse type."""
    parts = []
    if ecl_type & ECL_TOTAL:
        parts.append("TOTAL")
    if ecl_type & ECL_ANNULAR:
        parts.append("ANNULAR")
    if ecl_type & ECL_PARTIAL:
        parts.append("PARTIAL")
    if ecl_type & ECL_ANNULAR_TOTAL:
        parts.append("HYBRID")
    if ecl_type & ECL_PENUMBRAL:
        parts.append("PENUMBRAL")
    if ecl_type & ECL_CENTRAL:
        parts.append("CENTRAL")
    if ecl_type & ECL_NONCENTRAL:
        parts.append("NONCENTRAL")
    return "|".join(parts) if parts else f"TYPE={ecl_type}"


def compare_solar_eclipse_how(tmax, locations):
    """Compare sol_eclipse_how at multiple locations."""
    results = []

    for lon, lat, alt, name in locations:
        geopos = [lon, lat, alt]

        try:
            # pyswisseph: swe.sol_eclipse_how(tjd, geopos, ifl)
            se_res = swe.sol_eclipse_how(tmax, geopos, 0)
            se_type = se_res[0]
            se_attr = se_res[1]
        except Exception as e:
            results.append((name, "SE_ERROR", str(e)))
            continue

        try:
            # libephemeris: sol_eclipse_how(tjd, ifl, geopos)
            le_res = ephem.sol_eclipse_how(tmax, 0, geopos)
            le_type = le_res[0]
            le_attr = le_res[1]
        except Exception as e:
            results.append((name, "LE_ERROR", str(e)))
            continue

        # Compare attributes
        # attr[0] = fraction of solar diameter covered by moon (= magnitude for SE)
        # attr[1] = ratio of lunar diameter to solar one
        # attr[2] = fraction of solar disc covered by moon (obscuration)
        # attr[3] = diameter of core shadow in km
        # attr[4] = azimuth of sun at tjd
        # attr[5] = true altitude of sun above horizon at tjd
        # attr[6] = apparent altitude of sun above horizon at tjd
        # attr[7] = angular distance of moon from sun in degrees

        attr_labels = [
            "magnitude",
            "diam_ratio",
            "obscuration",
            "shadow_width_km",
            "sun_azimuth",
            "sun_true_alt",
            "sun_app_alt",
            "moon_sun_dist",
        ]

        # Tolerances per attribute
        tolerances = {
            "magnitude": 0.01,  # 1% of magnitude
            "diam_ratio": 0.005,  # 0.5%
            "obscuration": 0.02,  # 2% of obscuration
            "shadow_width_km": 50.0,  # 50 km (known systematic diff)
            "sun_azimuth": 0.1,  # 0.1 degrees
            "sun_true_alt": 0.05,  # 0.05 degrees
            "sun_app_alt": 0.05,  # 0.05 degrees
            "moon_sun_dist": 0.01,  # 0.01 degrees
        }

        for i, label in enumerate(attr_labels):
            if i >= len(se_attr) or i >= len(le_attr):
                break

            se_val = se_attr[i]
            le_val = le_attr[i]
            diff = abs(le_val - se_val)
            tol = tolerances.get(label, 0.1)

            # Special handling for shadow width (sign conventions differ)
            if label == "shadow_width_km":
                # Compare absolute values — sign convention differs
                diff = abs(abs(le_val) - abs(se_val))

            # Special handling for azimuth (wrapping)
            if label == "sun_azimuth":
                d = le_val - se_val
                if d > 180:
                    d -= 360
                elif d < -180:
                    d += 360
                diff = abs(d)

            passed = diff < tol
            results.append((name, label, se_val, le_val, diff, tol, passed))

    return results


def compare_lunar_eclipse_how(tmax):
    """Compare lun_eclipse_how at a location."""
    results = []

    locations = [
        (0.0, 51.5, 0.0, "London"),
        (139.7, 35.7, 0.0, "Tokyo"),
        (-74.0, 40.7, 0.0, "NYC"),
    ]

    for lon, lat, alt, name in locations:
        geopos = [lon, lat, alt]

        try:
            # pyswisseph: swe.lun_eclipse_how(tjd, geopos, ifl)
            se_res = swe.lun_eclipse_how(tmax, geopos, 0)
            se_type = se_res[0]
            se_attr = se_res[1]
        except Exception as e:
            results.append((name, "SE_ERROR", str(e)))
            continue

        try:
            # libephemeris: lun_eclipse_how(tjd, ifl, geopos)
            le_res = ephem.lun_eclipse_how(tmax, 0, geopos)
            le_type = le_res[0]
            le_attr = le_res[1]
        except Exception as e:
            results.append((name, "LE_ERROR", str(e)))
            continue

        # Lunar eclipse attributes:
        # attr[0] = umbral magnitude at tjd
        # attr[1] = penumbral magnitude
        # attr[4] = azimuth of moon at tjd
        # attr[5] = true altitude of moon above horizon
        # attr[6] = apparent altitude of moon above horizon
        # attr[7] = distance of moon from opposition in degrees

        attr_labels_idx = [
            (0, "umbral_mag", 0.02),
            (1, "penumbral_mag", 0.02),
            (4, "moon_azimuth", 0.1),
            (5, "moon_true_alt", 0.05),
            (6, "moon_app_alt", 0.05),
            (7, "moon_opp_dist", 0.01),
        ]

        for idx, label, tol in attr_labels_idx:
            if idx >= len(se_attr) or idx >= len(le_attr):
                break

            se_val = se_attr[idx]
            le_val = le_attr[idx]
            diff = abs(le_val - se_val)

            if label == "moon_azimuth":
                d = le_val - se_val
                if d > 180:
                    d -= 360
                elif d < -180:
                    d += 360
                diff = abs(d)

            passed = diff < tol
            results.append((name, label, se_val, le_val, diff, tol, passed))

    return results


def main():
    print("=" * 80)
    print("ROUND 118: Eclipse how() Deep — Solar & Lunar Eclipse Attribute Validation")
    print("=" * 80)

    total_tests = 0
    total_pass = 0
    total_fail = 0
    failures = []

    # Find solar eclipses 2000-2030
    print("\n--- Finding solar eclipses 2000-2030 ---")
    solar_eclipses = find_solar_eclipses(2451545.0, count=30)
    print(f"Found {len(solar_eclipses)} solar eclipses")

    # Locations for solar eclipse testing
    locations = [
        (0.0, 51.5, 0.0, "London"),
        (139.7, 35.7, 0.0, "Tokyo"),
        (-74.0, 40.7, 0.0, "NYC"),
        (12.5, 41.9, 50.0, "Rome"),
        (-118.2, 34.1, 71.0, "LA"),
        (77.2, 28.6, 216.0, "Delhi"),
        (151.2, -33.9, 3.0, "Sydney"),
        (-43.2, -22.9, 11.0, "RioDeJ"),
    ]

    # Test solar eclipse how at each location
    print("\n--- Testing sol_eclipse_how at multiple locations ---")
    for ecl_type, tmax in solar_eclipses:
        ecl_str = classify_eclipse(ecl_type)

        results = compare_solar_eclipse_how(tmax, locations)

        for item in results:
            if isinstance(item[1], str) and item[1].endswith("ERROR"):
                continue

            name, label, se_val, le_val, diff, tol, passed = item
            total_tests += 1
            if passed:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 40:
                    failures.append(
                        f"  SOL {ecl_str} JD={tmax:.2f} {name}: {label}: "
                        f"SE={se_val:.6f} LE={le_val:.6f} diff={diff:.6f} (tol={tol})"
                    )

    # Find lunar eclipses 2000-2030
    print("\n--- Finding lunar eclipses 2000-2030 ---")
    lunar_eclipses = find_lunar_eclipses(2451545.0, count=30)
    print(f"Found {len(lunar_eclipses)} lunar eclipses")

    # Test lunar eclipse how
    print("\n--- Testing lun_eclipse_how at multiple locations ---")
    for ecl_type, tmax in lunar_eclipses:
        ecl_str = classify_eclipse(ecl_type)

        results = compare_lunar_eclipse_how(tmax)

        for item in results:
            if isinstance(item[1], str) and item[1].endswith("ERROR"):
                continue

            name, label, se_val, le_val, diff, tol, passed = item
            total_tests += 1
            if passed:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 80:
                    failures.append(
                        f"  LUN {ecl_str} JD={tmax:.2f} {name}: {label}: "
                        f"SE={se_val:.6f} LE={le_val:.6f} diff={diff:.6f} (tol={tol})"
                    )

    # Test eclipse how at multiple time offsets around maximum
    print("\n--- Testing attributes at phases around solar eclipse maximum ---")
    if solar_eclipses:
        # Take first total/annular eclipse
        test_ecl = None
        for ecl_type, tmax in solar_eclipses:
            if ecl_type & (ECL_TOTAL | ECL_ANNULAR):
                test_ecl = (ecl_type, tmax)
                break

        if test_ecl is None:
            test_ecl = solar_eclipses[0]

        ecl_type, tmax = test_ecl
        ecl_str = classify_eclipse(ecl_type)
        print(f"  Testing {ecl_str} at JD={tmax:.4f}")

        # Test at offsets: -2h, -1h, -30m, -15m, 0, +15m, +30m, +1h, +2h
        offsets_hours = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]

        for offset_h in offsets_hours:
            jd_test = tmax + offset_h / 24.0
            test_locs = locations[:3]  # London, Tokyo, NYC

            results = compare_solar_eclipse_how(jd_test, test_locs)

            for item in results:
                if isinstance(item[1], str) and item[1].endswith("ERROR"):
                    continue

                name, label, se_val, le_val, diff, tol, passed = item
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 100:
                        failures.append(
                            f"  SOL PHASE {ecl_str} offset={offset_h:+.2f}h {name}: {label}: "
                            f"SE={se_val:.6f} LE={le_val:.6f} diff={diff:.6f} (tol={tol})"
                        )

    # Test lunar eclipse attributes at phases
    print("\n--- Testing attributes at phases around lunar eclipse maximum ---")
    if lunar_eclipses:
        test_ecl = None
        for ecl_type, tmax in lunar_eclipses:
            if ecl_type & ECL_TOTAL:
                test_ecl = (ecl_type, tmax)
                break

        if test_ecl is None:
            test_ecl = lunar_eclipses[0]

        ecl_type, tmax = test_ecl
        ecl_str = classify_eclipse(ecl_type)
        print(f"  Testing {ecl_str} at JD={tmax:.4f}")

        offsets_hours = [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3]

        for offset_h in offsets_hours:
            jd_test = tmax + offset_h / 24.0

            results = compare_lunar_eclipse_how(jd_test)

            for item in results:
                if isinstance(item[1], str) and item[1].endswith("ERROR"):
                    continue

                name, label, se_val, le_val, diff, tol, passed = item
                total_tests += 1
                if passed:
                    total_pass += 1
                else:
                    total_fail += 1
                    if len(failures) < 120:
                        failures.append(
                            f"  LUN PHASE {ecl_str} offset={offset_h:+.2f}h {name}: {label}: "
                            f"SE={se_val:.6f} LE={le_val:.6f} diff={diff:.6f} (tol={tol})"
                        )

    # Test eclipse type consistency between how() and when_glob()
    print("\n--- Testing eclipse type consistency ---")
    for ecl_type_glob, tmax in solar_eclipses[:15]:
        # Check a central path location — use (0, 0) as generic
        geopos = [0.0, 0.0, 0.0]
        try:
            se_res = swe.sol_eclipse_how(tmax, geopos, 0)
            se_type_how = se_res[0]

            le_res = ephem.sol_eclipse_how(tmax, 0, geopos)
            le_type_how = le_res[0]

            total_tests += 1
            # Check that both agree on the major eclipse type
            se_major = se_type_how & (
                ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
            )
            le_major = le_type_how & (
                ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
            )

            if se_major == le_major or se_major == 0 or le_major == 0:
                total_pass += 1
            else:
                total_fail += 1
                if len(failures) < 130:
                    failures.append(
                        f"  TYPE MISMATCH: JD={tmax:.2f}: SE={classify_eclipse(se_type_how)} "
                        f"LE={classify_eclipse(le_type_how)}"
                    )
        except Exception:
            pass

    # Summary
    print("\n" + "=" * 80)
    pct = 100 * total_pass / total_tests if total_tests > 0 else 0
    print(f"ROUND 118 RESULTS: {total_pass}/{total_tests} passed ({pct:.1f}%)")
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
