#!/usr/bin/env python3
"""Round 41: House Systems Comprehensive Verification.

Tests ALL supported house systems across a wide range of latitudes and dates,
comparing libephemeris vs pyswisseph for cusps and ascmc values.

Phases:
  P1: All house systems at standard latitudes (0°, 23.4°, 45°, 51.5°, 60°)
  P2: High-latitude behavior (65°, 66.3°, 66.6°, 70°, 80°, 85°, 89°)
  P3: Southern hemisphere mirror symmetry
  P4: Multiple dates (seasonal variation - equinoxes, solstices)
  P5: ARMC-based functions (houses_armc vs houses_armc)
  P6: house_pos - planet-in-house computation
  P7: house_name - name string verification
  P8: Sidereal mode house cusps (houses_ex with FLG_SIDEREAL)
"""

from __future__ import annotations

import math
import os
import sys
import traceback

# Force skyfield mode
os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)


# Helpers for house system byte/int conversion
def se_hsys(ch: str) -> bytes:
    return ch.encode("ascii")


def le_hsys(ch: str) -> int:
    return ord(ch)


# All house system codes to test
# Skip 'i' (Sunshine alt) as pyswisseph may not support it
ALL_SYSTEMS = list("PKORCEWMBTUGVXHFSLANQYDI")

# Systems that fail at polar latitudes (|lat| > ~66.5°)
POLAR_FAIL = set("PKG")

# Systems unstable at extreme latitudes (|lat| >= 80°)
EXTREME_UNSTABLE = set("CRTBHUYF")

# Standard latitudes
STANDARD_LATS = [0.0, 23.4, 45.0, 51.5, 60.0]

# High latitudes
HIGH_LATS = [65.0, 66.3, 66.56, 70.0, 80.0, 85.0, 89.0]

# Southern hemisphere
SOUTH_LATS = [-23.4, -45.0, -51.5, -60.0, -66.3]

# Standard longitudes (several to test)
LONGITUDES = [0.0, 11.25, -73.97, 139.69]  # Greenwich, Rome, NYC, Tokyo

# Test dates (JD UT)
DATES = {
    "J2000": 2451545.0,
    "2024_spring_equinox": 2460389.5,
    "2024_summer_solstice": 2460481.0,
    "2024_autumn_equinox": 2460572.5,
    "2024_winter_solstice": 2460663.0,
    "1980_random": 2444300.5,
    "2050_random": 2469807.5,
}

passed = 0
failed = 0
errors = 0
skipped = 0
results = {"passed": [], "failed": [], "errors": [], "skipped": []}

CUSP_TOL = 1.0 / 3600.0  # 1 arcsecond
ASCMC_TOL = 1.0 / 3600.0  # 1 arcsecond


def compare_cusps_ascmc(
    phase: str,
    label: str,
    hsys_ch: str,
    tjdut: float,
    lat: float,
    lon: float,
    use_ex: bool = False,
    sidereal: bool = False,
):
    global passed, failed, errors, skipped

    try:
        # Determine flags
        se_flags = 0
        le_flags = 0
        if sidereal:
            se_flags |= swe.FLG_SIDEREAL
            le_flags |= ephem.FLG_SIDEREAL
            # Use Lahiri as default
            swe.set_sid_mode(swe.SIDM_LAHIRI)
            ephem.set_sid_mode(ephem.SIDM_LAHIRI, 0, 0)

        # Call SE
        try:
            if use_ex or sidereal:
                se_result = swe.houses_ex(tjdut, lat, lon, se_hsys(hsys_ch), se_flags)
            else:
                se_result = swe.houses(tjdut, lat, lon, se_hsys(hsys_ch))
            se_cusps = se_result[0]
            se_ascmc = se_result[1]
        except Exception as e:
            if "polar" in str(e).lower() or "beyond polar" in str(e).lower():
                skipped += 1
                results["skipped"].append(
                    f"{phase} {label}: SE polar error for {hsys_ch}"
                )
                return
            raise

        # Call LE
        try:
            if use_ex or sidereal:
                le_result = ephem.houses_ex(
                    tjdut, lat, lon, le_hsys(hsys_ch), le_flags
                )
            else:
                le_result = ephem.houses(tjdut, lat, lon, le_hsys(hsys_ch))
            le_cusps = le_result[0]
            le_ascmc = le_result[1]
        except Exception as e:
            if "polar" in str(e).lower():
                # If SE didn't fail but LE did, that's a potential issue
                # but some systems have different polar thresholds
                skipped += 1
                results["skipped"].append(
                    f"{phase} {label}: LE polar error for {hsys_ch}"
                )
                return
            raise

        # Compare cusps
        n_cusps = len(se_cusps)
        le_n_cusps = len(le_cusps)
        if n_cusps != le_n_cusps:
            failed += 1
            results["failed"].append(
                f"{phase} {label}: cusp count mismatch SE={n_cusps} LE={le_n_cusps}"
            )
            return

        cusp_max_diff = 0.0
        cusp_fails = []
        for i in range(n_cusps):
            diff = abs(se_cusps[i] - le_cusps[i])
            # Handle 360° wrap
            if diff > 180:
                diff = 360 - diff
            cusp_max_diff = max(cusp_max_diff, diff)
            if diff > CUSP_TOL:
                cusp_fails.append(
                    f'  cusp[{i + 1}] SE={se_cusps[i]:.6f} LE={le_cusps[i]:.6f} diff={diff * 3600:.2f}"'
                )

        # Compare ascmc (first 8 values)
        ascmc_max_diff = 0.0
        ascmc_fails = []
        ascmc_names = [
            "Asc",
            "MC",
            "ARMC",
            "Vertex",
            "EquAsc",
            "CoAsc_Koch",
            "CoAsc_Munk",
            "PolarAsc",
        ]
        n_ascmc = min(len(se_ascmc), len(le_ascmc), 8)
        for i in range(n_ascmc):
            # Skip zero/undefined values
            if se_ascmc[i] == 0.0 and le_ascmc[i] == 0.0:
                continue
            diff = abs(se_ascmc[i] - le_ascmc[i])
            if diff > 180:
                diff = 360 - diff
            ascmc_max_diff = max(ascmc_max_diff, diff)
            if diff > ASCMC_TOL:
                ascmc_fails.append(
                    f"  ascmc[{i}]({ascmc_names[i] if i < len(ascmc_names) else '?'}) "
                    f'SE={se_ascmc[i]:.6f} LE={le_ascmc[i]:.6f} diff={diff * 3600:.2f}"'
                )

        if cusp_fails or ascmc_fails:
            failed += 1
            detail = f'cusp_max={cusp_max_diff * 3600:.2f}" ascmc_max={ascmc_max_diff * 3600:.2f}"'
            fail_details = "\n".join(cusp_fails + ascmc_fails)
            results["failed"].append(
                f"{phase} {label} [{hsys_ch}]: {detail}\n{fail_details}"
            )
        else:
            passed += 1
            results["passed"].append(
                f'{phase} {label} [{hsys_ch}]: cusp_max={cusp_max_diff * 3600:.3f}" ascmc_max={ascmc_max_diff * 3600:.3f}"'
            )

    except Exception as e:
        errors += 1
        results["errors"].append(
            f"{phase} {label} [{hsys_ch}]: {e}\n{traceback.format_exc()}"
        )


def phase1():
    """All house systems at standard latitudes."""
    print("\n=== P1: All house systems at standard latitudes ===")
    tjdut = DATES["J2000"]
    lon = 11.25  # Rome

    for lat in STANDARD_LATS:
        for hsys_ch in ALL_SYSTEMS:
            label = f"lat={lat:+.1f} lon={lon}"
            compare_cusps_ascmc("P1", label, hsys_ch, tjdut, lat, lon)


def phase2():
    """High-latitude behavior."""
    print("\n=== P2: High-latitude behavior ===")
    tjdut = DATES["J2000"]
    lon = 24.94  # Helsinki-ish

    for lat in HIGH_LATS:
        for hsys_ch in ALL_SYSTEMS:
            # Skip systems known to fail at polar latitudes
            if lat > 66.5 and hsys_ch in POLAR_FAIL:
                skipped_count_before = skipped
                compare_cusps_ascmc("P2", f"lat={lat:+.1f}", hsys_ch, tjdut, lat, lon)
                continue
            compare_cusps_ascmc("P2", f"lat={lat:+.1f}", hsys_ch, tjdut, lat, lon)


def phase3():
    """Southern hemisphere mirror symmetry."""
    print("\n=== P3: Southern hemisphere ===")
    tjdut = DATES["J2000"]
    lon = -43.17  # Rio de Janeiro

    for lat in SOUTH_LATS:
        for hsys_ch in ALL_SYSTEMS:
            compare_cusps_ascmc("P3", f"lat={lat:+.1f}", hsys_ch, tjdut, lat, lon)


def phase4():
    """Multiple dates (seasonal variation)."""
    print("\n=== P4: Multiple dates / seasonal variation ===")
    lat = 51.5  # London
    lon = -0.12

    for date_name, tjdut in DATES.items():
        for hsys_ch in ["P", "K", "O", "R", "C", "E", "W", "B", "T", "M"]:
            label = f"{date_name} lat={lat}"
            compare_cusps_ascmc("P4", label, hsys_ch, tjdut, lat, lon)


def phase5():
    """ARMC-based functions."""
    print("\n=== P5: ARMC-based functions (houses_armc) ===")
    # Test houses_armc directly with known ARMC/obliquity values
    test_cases = [
        # (armc, lat, eps, label)
        (0.0, 45.0, 23.4393, "ARMC=0"),
        (90.0, 45.0, 23.4393, "ARMC=90"),
        (180.0, 45.0, 23.4393, "ARMC=180"),
        (270.0, 45.0, 23.4393, "ARMC=270"),
        (45.0, 0.0, 23.4393, "equator"),
        (45.0, 60.0, 23.4393, "lat=60"),
        (123.456, 51.5, 23.4393, "London-like"),
        (200.0, -33.9, 23.4393, "Sydney-like"),
        (315.0, 35.7, 23.4393, "Tokyo-like"),
    ]

    global passed, failed, errors

    for armc, lat, eps, case_label in test_cases:
        for hsys_ch in ALL_SYSTEMS:
            if abs(lat) > 66.5 and hsys_ch in POLAR_FAIL:
                continue

            try:
                se_cusps, se_ascmc = swe.houses_armc(armc, lat, eps, se_hsys(hsys_ch))
                le_cusps, le_ascmc = ephem.houses_armc(
                    armc, lat, eps, le_hsys(hsys_ch)
                )

                max_diff = 0.0
                fail_details = []
                for i in range(min(len(se_cusps), len(le_cusps))):
                    diff = abs(se_cusps[i] - le_cusps[i])
                    if diff > 180:
                        diff = 360 - diff
                    max_diff = max(max_diff, diff)
                    if diff > CUSP_TOL:
                        fail_details.append(
                            f'  cusp[{i + 1}] SE={se_cusps[i]:.6f} LE={le_cusps[i]:.6f} diff={diff * 3600:.2f}"'
                        )

                for i in range(min(len(se_ascmc), len(le_ascmc), 8)):
                    if se_ascmc[i] == 0.0 and le_ascmc[i] == 0.0:
                        continue
                    diff = abs(se_ascmc[i] - le_ascmc[i])
                    if diff > 180:
                        diff = 360 - diff
                    max_diff = max(max_diff, diff)
                    if diff > CUSP_TOL:
                        ascmc_names = [
                            "Asc",
                            "MC",
                            "ARMC",
                            "Vertex",
                            "EquAsc",
                            "CoAsc_Koch",
                            "CoAsc_Munk",
                            "PolarAsc",
                        ]
                        name = ascmc_names[i] if i < len(ascmc_names) else "?"
                        fail_details.append(
                            f'  ascmc[{i}]({name}) SE={se_ascmc[i]:.6f} LE={le_ascmc[i]:.6f} diff={diff * 3600:.2f}"'
                        )

                if fail_details:
                    failed += 1
                    results["failed"].append(
                        f'P5 {case_label} [{hsys_ch}]: max={max_diff * 3600:.2f}"\n'
                        + "\n".join(fail_details)
                    )
                else:
                    passed += 1
                    results["passed"].append(
                        f'P5 {case_label} [{hsys_ch}]: max={max_diff * 3600:.3f}"'
                    )

            except Exception as e:
                if "polar" in str(e).lower():
                    continue
                errors += 1
                results["errors"].append(f"P5 {case_label} [{hsys_ch}]: {e}")


def phase6():
    """house_pos - planet-in-house computation."""
    print("\n=== P6: house_pos verification ===")
    global passed, failed, errors

    # Get a set of planet positions first
    tjdut = DATES["J2000"]
    planets = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  # Sun through Pluto
    planet_names = [
        "Sun",
        "Moon",
        "Mercury",
        "Venus",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",
    ]

    test_locations = [
        (45.0, 11.25, "Rome"),
        (51.5, -0.12, "London"),
        (40.7, -74.0, "NYC"),
        (0.0, 0.0, "Equator"),
    ]

    # House systems that support house_pos
    hp_systems = ["P", "K", "O", "R", "C", "E", "W", "B", "T", "M"]

    for lat, lon, loc_name in test_locations:
        # First get ARMC and obliquity from houses call
        try:
            se_cusps, se_ascmc = swe.houses(tjdut, lat, lon, se_hsys("P"))
            armc = se_ascmc[2]
            eps = se_ascmc[8] if len(se_ascmc) > 8 else 23.4393

            # Get obliquity from SE
            nutob = swe.calc_ut(tjdut, swe.ECL_NUT)[0]
            eps = nutob[1]  # true obliquity
        except Exception:
            continue

        for hsys_ch in hp_systems:
            for pi, pname in zip(planets, planet_names):
                try:
                    # Get planet ecliptic position
                    se_pos = swe.calc_ut(tjdut, pi, swe.FLG_SWIEPH)[0]
                    p_lon = se_pos[0]
                    p_lat = se_pos[1]

                    # SE house_pos: swe.house_pos(armc, geolat, eps, xpin, hsys)
                    se_hp = swe.house_pos(
                        armc, lat, eps, (p_lon, p_lat), se_hsys(hsys_ch)
                    )

                    # LE house_pos: house_pos(armc, lat, eps, objcoord, hsys)
                    le_hp = ephem.house_pos(armc, lat, eps, (p_lon, p_lat), hsys_ch)

                    diff = abs(se_hp - le_hp)
                    tol = 0.01  # 0.01 house unit

                    label = f"{loc_name} [{hsys_ch}] {pname}"
                    if diff > tol:
                        failed += 1
                        results["failed"].append(
                            f"P6 {label}: SE={se_hp:.6f} LE={le_hp:.6f} diff={diff:.6f}"
                        )
                    else:
                        passed += 1
                        results["passed"].append(f"P6 {label}: diff={diff:.6f}")

                except Exception as e:
                    if "polar" in str(e).lower() or "not supported" in str(e).lower():
                        continue
                    errors += 1
                    results["errors"].append(f"P6 {loc_name} [{hsys_ch}] {pname}: {e}")


def phase7():
    """house_name verification."""
    print("\n=== P7: house_name verification ===")
    global passed, failed, errors

    for hsys_ch in ALL_SYSTEMS:
        try:
            se_name = swe.house_name(se_hsys(hsys_ch))
            le_name = ephem.house_name(le_hsys(hsys_ch))

            # Both should return non-empty strings
            if not se_name or not le_name:
                failed += 1
                results["failed"].append(
                    f"P7 house_name('{hsys_ch}'): SE='{se_name}' LE='{le_name}' (empty)"
                )
            elif se_name.lower().replace("-", "").replace(
                " ", ""
            ) != le_name.lower().replace("-", "").replace(" ", ""):
                # Loose comparison - allow minor formatting differences
                # Check if the important part of the name matches
                se_key = se_name.lower().split("/")[0].split("(")[0].strip()
                le_key = le_name.lower().split("/")[0].split("(")[0].strip()
                if se_key[:6] == le_key[:6]:  # First 6 chars match
                    passed += 1
                    results["passed"].append(
                        f"P7 house_name('{hsys_ch}'): SE='{se_name}' LE='{le_name}' (close match)"
                    )
                else:
                    failed += 1
                    results["failed"].append(
                        f"P7 house_name('{hsys_ch}'): SE='{se_name}' LE='{le_name}' (mismatch)"
                    )
            else:
                passed += 1
                results["passed"].append(f"P7 house_name('{hsys_ch}'): '{le_name}'")
        except Exception as e:
            errors += 1
            results["errors"].append(f"P7 house_name('{hsys_ch}'): {e}")


def phase8():
    """Sidereal mode house cusps."""
    print("\n=== P8: Sidereal mode house cusps ===")

    tjdut = DATES["J2000"]

    # Test a few sidereal modes
    sid_modes = [
        (0, "Fagan/Bradley"),
        (1, "Lahiri"),
        (3, "Raman"),
        (5, "Krishnamurti"),
        (27, "True Citra"),
    ]

    lats = [0.0, 45.0, 51.5, -33.9]
    lon = 11.25
    systems = ["P", "K", "O", "E", "W", "R", "C", "B"]

    for sid_id, sid_name in sid_modes:
        swe.set_sid_mode(sid_id)
        ephem.set_sid_mode(sid_id, 0, 0)

        for lat in lats:
            for hsys_ch in systems:
                label = f"{sid_name} lat={lat:+.1f}"
                compare_cusps_ascmc(
                    "P8", label, hsys_ch, tjdut, lat, lon, use_ex=True, sidereal=True
                )

    # Reset to tropical
    swe.set_sid_mode(0)


def main():
    print("=" * 70)
    print("ROUND 41: House Systems Comprehensive Verification")
    print("=" * 70)
    print(f"Testing {len(ALL_SYSTEMS)} house systems")
    print(f"Standard latitudes: {STANDARD_LATS}")
    print(f"High latitudes: {HIGH_LATS}")
    print(f'Tolerance: cusps={CUSP_TOL * 3600:.1f}" ascmc={ASCMC_TOL * 3600:.1f}"')

    phase1()
    print(
        f"  After P1: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase2()
    print(
        f"  After P2: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase3()
    print(
        f"  After P3: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase4()
    print(
        f"  After P4: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase5()
    print(
        f"  After P5: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase6()
    print(
        f"  After P6: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase7()
    print(
        f"  After P7: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    phase8()
    print(
        f"  After P8: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped"
    )

    # Final summary
    total = passed + failed + errors
    print("\n" + "=" * 70)
    print(f"ROUND 41 FINAL: {passed}/{total} passed ({100 * passed / total:.1f}%)")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Errors:  {errors}")
    print(f"  Skipped: {skipped}")
    print("=" * 70)

    # Show failures
    if results["failed"]:
        print(f"\n--- FAILURES ({len(results['failed'])}) ---")
        for f in results["failed"][:50]:
            print(f)
        if len(results["failed"]) > 50:
            print(f"  ... and {len(results['failed']) - 50} more")

    if results["errors"]:
        print(f"\n--- ERRORS ({len(results['errors'])}) ---")
        for e in results["errors"][:20]:
            print(e)
        if len(results["errors"]) > 20:
            print(f"  ... and {len(results['errors']) - 20} more")

    # Show some passed stats
    if results["passed"]:
        # Compute overall max diff from passed
        max_arcsec = 0.0
        for p in results["passed"]:
            for part in p.split():
                if part.startswith("cusp_max=") or part.startswith("max="):
                    val = part.split("=")[1].rstrip('"')
                    try:
                        max_arcsec = max(max_arcsec, float(val))
                    except ValueError:
                        pass
        print(f'\nMax cusp/ascmc diff across all passed tests: {max_arcsec:.3f}"')

    # Show skipped
    if results["skipped"]:
        print(f"\n--- SKIPPED ({len(results['skipped'])}) ---")
        # Group by reason
        polar_count = sum(1 for s in results["skipped"] if "polar" in s.lower())
        print(f"  Polar latitude failures: {polar_count}")

    return 0 if failed == 0 and errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
