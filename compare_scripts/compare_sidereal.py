"""
Sidereal/Tropical Modes Comparison Script

Compares all ayanamsha (sidereal) modes between pyswisseph and libephemeris.
Tests all 43 ayanamsha systems across different dates and planets.
"""

import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
from comparison_utils import (
    angular_diff,
    format_coord,
    format_diff,
    format_status,
    TestStatistics,
    print_header,
    parse_args,
    STANDARD_SUBJECTS,
    Tolerances,
)
import sys

# ============================================================================
# AYANAMSHA MODES TO TEST
# ============================================================================

AYANAMSHA_MODES = {
    SE_SIDM_FAGAN_BRADLEY: "Fagan/Bradley",
    SE_SIDM_LAHIRI: "Lahiri",
    SE_SIDM_DELUCE: "De Luce",
    SE_SIDM_RAMAN: "Raman",
    SE_SIDM_USHASHASHI: "Ushashashi",
    SE_SIDM_KRISHNAMURTI: "Krishnamurti",
    SE_SIDM_DJWHAL_KHUL: "Djwhal Khul",
    SE_SIDM_YUKTESHWAR: "Yukteshwar",
    SE_SIDM_JN_BHASIN: "JN Bhasin",
    SE_SIDM_BABYL_KUGLER1: "Babylonian/Kugler 1",
    SE_SIDM_BABYL_KUGLER2: "Babylonian/Kugler 2",
    SE_SIDM_BABYL_KUGLER3: "Babylonian/Kugler 3",
    SE_SIDM_BABYL_HUBER: "Babylonian/Huber",
    SE_SIDM_BABYL_ETPSC: "Babylonian/ETPSC",
    SE_SIDM_ALDEBARAN_15TAU: "Aldebaran at 15 Tau",
    SE_SIDM_HIPPARCHOS: "Hipparchos",
    SE_SIDM_SASSANIAN: "Sassanian",
    SE_SIDM_GALCENT_0SAG: "Galactic Center at 0 Sag",
    SE_SIDM_J2000: "J2000",
    SE_SIDM_J1900: "J1900",
    SE_SIDM_B1950: "B1950",
    SE_SIDM_SURYASIDDHANTA: "Suryasiddhanta",
    SE_SIDM_SURYASIDDHANTA_MSUN: "Suryasiddhanta (mean Sun)",
    SE_SIDM_ARYABHATA: "Aryabhata",
    SE_SIDM_ARYABHATA_MSUN: "Aryabhata (mean Sun)",
    SE_SIDM_SS_REVATI: "SS Revati",
    SE_SIDM_SS_CITRA: "SS Citra",
    SE_SIDM_TRUE_CITRA: "True Citra",
    SE_SIDM_TRUE_REVATI: "True Revati",
    SE_SIDM_TRUE_PUSHYA: "True Pushya",
    SE_SIDM_GALCENT_RGILBRAND: "Galactic Center (Gil Brand)",
    SE_SIDM_GALEQU_IAU1958: "Galactic Equator (IAU 1958)",
    SE_SIDM_GALEQU_TRUE: "Galactic Equator (True)",
    SE_SIDM_GALEQU_MULA: "Galactic Equator at Mula",
    SE_SIDM_GALALIGN_MARDYKS: "Galactic Alignment (Mardyks)",
    SE_SIDM_TRUE_MULA: "True Mula",
    SE_SIDM_GALCENT_MULA_WILHELM: "Galactic Center at Mula (Wilhelm)",
    SE_SIDM_ARYABHATA_522: "Aryabhata 522",
    SE_SIDM_BABYL_BRITTON: "Babylonian (Britton)",
    SE_SIDM_TRUE_SHEORAN: "True Sheoran",
    SE_SIDM_GALCENT_COCHRANE: "Galactic Center (Cochrane)",
    SE_SIDM_GALEQU_FIORENZA: "Galactic Equator (Fiorenza)",
    SE_SIDM_VALENS_MOON: "Valens (Moon)",
}

# Test planets (fast and slow movers)
TEST_PLANETS = {
    SE_SUN: "Sun",
    SE_MOON: "Moon",
    SE_MARS: "Mars",
    SE_JUPITER: "Jupiter",
    SE_SATURN: "Saturn",
}

# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_sidereal_position(
    subject_name: str,
    date_str: str,
    jd: float,
    planet_id: int,
    planet_name: str,
    sid_mode: int,
    sid_name: str,
    verbose: bool = False,
) -> tuple:
    """
    Compare sidereal position calculation for a specific ayanamsha mode.

    Returns:
        (passed, max_diff, error_occurred)
    """
    # Set sidereal mode
    swe.set_sid_mode(sid_mode)
    pyephem.swe_set_sid_mode(sid_mode)

    flags = SEFLG_SIDEREAL | SEFLG_SPEED

    # Calculate with SwissEph
    try:
        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        lon_swe, lat_swe, dist_swe = res_swe[0], res_swe[1], res_swe[2]
        lon_speed_swe, lat_speed_swe, dist_speed_swe = (
            res_swe[3],
            res_swe[4],
            res_swe[5],
        )
    except Exception as e:
        if verbose:
            print(
                f"[{subject_name}] [{date_str}] [{sid_name:<35}] {planet_name}: SWE ERROR {e}"
            )
        return False, 0.0, True

    # Calculate with Python Ephemeris
    try:
        res_py, _ = pyephem.swe_calc_ut(jd, planet_id, flags)
        lon_py, lat_py, dist_py = res_py[0], res_py[1], res_py[2]
        lon_speed_py, lat_speed_py, dist_speed_py = res_py[3], res_py[4], res_py[5]
    except Exception as e:
        if verbose:
            print(
                f"[{subject_name}] [{date_str}] [{sid_name:<35}] {planet_name}: PY ERROR {e}"
            )
        return False, 0.0, True

    # Calculate differences
    diff_lon = angular_diff(lon_swe, lon_py)
    diff_lat = abs(lat_swe - lat_py)
    diff_dist = abs(dist_swe - dist_py)
    diff_lon_speed = abs(lon_speed_swe - lon_speed_py)
    diff_lat_speed = abs(lat_speed_swe - lat_speed_py)
    diff_dist_speed = abs(dist_speed_swe - dist_speed_py)

    max_diff = max(diff_lon, diff_lat)

    # Check tolerances (use relaxed for position, strict for velocity)
    passed = (
        diff_lon < Tolerances.LONGITUDE_STRICT
        and diff_lat < Tolerances.LATITUDE_STRICT
        and diff_dist < Tolerances.DISTANCE_STRICT
        and diff_lon_speed < Tolerances.VELOCITY_ANGULAR
        and diff_lat_speed < Tolerances.VELOCITY_ANGULAR
        and diff_dist_speed < Tolerances.VELOCITY_RADIAL
    )

    status = format_status(passed)

    if verbose:
        print(f"\n{'=' * 80}")
        print(f"{subject_name} - {date_str} - {sid_name} - {planet_name}")
        print(f"{'=' * 80}")
        print("\nPosition:")
        print(
            f"  Longitude:  SWE={format_coord(lon_swe)}°  PY={format_coord(lon_py)}°  Diff={format_diff(diff_lon, 6)}°"
        )
        print(
            f"  Latitude:   SWE={format_coord(lat_swe)}°  PY={format_coord(lat_py)}°  Diff={format_diff(diff_lat, 6)}°"
        )
        print(
            f"  Distance:   SWE={format_coord(dist_swe)} AU  PY={format_coord(dist_py)} AU  Diff={format_diff(diff_dist, 6)} AU"
        )
        print("\nVelocity:")
        print(
            f"  Lon Speed:  SWE={format_coord(lon_speed_swe, 8)}°/d  PY={format_coord(lon_speed_py, 8)}°/d  Diff={format_diff(diff_lon_speed, 6)}"
        )
        print(
            f"  Lat Speed:  SWE={format_coord(lat_speed_swe, 8)}°/d  PY={format_coord(lat_speed_py, 8)}°/d  Diff={format_diff(diff_lat_speed, 6)}"
        )
        print(
            f"  Dist Speed: SWE={format_coord(dist_speed_swe, 8)} AU/d  PY={format_coord(dist_speed_py, 8)} AU/d  Diff={format_diff(diff_dist_speed, 6)}"
        )
        print(f"\nStatus: {'PASSED ✓' if passed else 'FAILED ✗'}")
    else:
        # Single line with all position and velocity values
        print(
            f"[{subject_name}] [{date_str}] [{sid_name:<35}] {planet_name:<10} "
            f"Lon={format_coord(lon_swe, 4, 8)}/{format_coord(lon_py, 4, 8)} "
            f"Lat={format_coord(lat_swe, 4, 8)}/{format_coord(lat_py, 4, 8)} "
            f"Dist={format_coord(dist_swe, 4, 8)}/{format_coord(dist_py, 4, 8)} "
            f"DiffL={format_diff(diff_lon, 4, 6)} DiffB={format_diff(diff_lat, 4, 6)} "
            f"SpeedL={format_coord(lon_speed_swe, 4, 10)}/{format_coord(lon_speed_py, 4, 10)} {status}"
        )

    return passed, max_diff, False


# ============================================================================
# AYANAMSHA VALUE COMPARISON
# ============================================================================


def compare_ayanamsha_value(
    date_str: str, jd: float, sid_mode: int, sid_name: str, verbose: bool = False
) -> tuple:
    """
    Compare ayanamsha value for a specific mode.

    Returns:
        (passed, diff, error_occurred)
    """
    # Set sidereal mode
    swe.set_sid_mode(sid_mode)
    pyephem.swe_set_sid_mode(sid_mode)

    # Get ayanamsha values
    try:
        ayan_swe = swe.get_ayanamsa_ut(jd)
        ayan_py = pyephem.swe_get_ayanamsa_ut(jd)
    except Exception as e:
        if verbose:
            print(f"[{date_str}] [{sid_name:<35}]: ERROR {e}")
        return False, 0.0, True

    diff = abs(ayan_swe - ayan_py)
    passed = diff < Tolerances.AYANAMSHA  # Use ayanamsha-specific tolerance
    status = format_status(passed)

    if verbose:
        print(
            f"[{date_str}] [{sid_name:<35}] Ayanamsha: SWE={format_coord(ayan_swe)}° PY={format_coord(ayan_py)}° "
            f"Diff={format_diff(diff, 6)}° {status}"
        )
    else:
        # Single line with both values
        print(
            f"[{date_str}] [{sid_name:<35}] Ayan={format_coord(ayan_swe, 4, 8)}/{format_coord(ayan_py, 4, 8)} "
            f"Diff={format_diff(diff, 4, 6)} {status}"
        )

    return passed, diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(
    verbose: bool = False, planets_only: bool = False, ayanamsha_only: bool = False
) -> tuple:
    """
    Run all sidereal mode comparison tests.

    Args:
        verbose: If True, print detailed output
        planets_only: If True, only test planetary positions
        ayanamsha_only: If True, only test ayanamsha values

    Returns:
        (passed_count, total_count)
    """
    print_header("SIDEREAL/TROPICAL MODES COMPARISON")

    stats = TestStatistics()

    # Test ayanamsha values if requested
    if not planets_only:
        print("\n--- Ayanamsha Values ---\n")
        for sid_mode, sid_name in AYANAMSHA_MODES.items():
            for name, year, month, day, hour, lat, lon, alt in STANDARD_SUBJECTS[
                :3
            ]:  # Use first 3 subjects
                jd = swe.julday(year, month, day, hour)
                date_str = f"{year}-{month:02d}-{day:02d}"

                passed, diff, error = compare_ayanamsha_value(
                    date_str=date_str,
                    jd=jd,
                    sid_mode=sid_mode,
                    sid_name=sid_name,
                    verbose=verbose,
                )

                stats.add_result(passed, diff, error)

    # Test planetary positions if requested
    if not ayanamsha_only:
        print("\n--- Sidereal Planetary Positions ---\n")
        for name, year, month, day, hour, lat, lon, alt in STANDARD_SUBJECTS[
            :2
        ]:  # Use first 2 subjects
            jd = swe.julday(year, month, day, hour)
            date_str = f"{year}-{month:02d}-{day:02d}"

            for sid_mode, sid_name in list(AYANAMSHA_MODES.items())[
                :10
            ]:  # Test first 10 modes
                for planet_id, planet_name in TEST_PLANETS.items():
                    passed, diff, error = compare_sidereal_position(
                        subject_name=name,
                        date_str=date_str,
                        jd=jd,
                        planet_id=planet_id,
                        planet_name=planet_name,
                        sid_mode=sid_mode,
                        sid_name=sid_name,
                        verbose=verbose,
                    )

                    stats.add_result(passed, diff, error)

    # Print summary
    stats.print_summary("SIDEREAL MODES COMPARISON SUMMARY")

    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_sidereal.py [OPTIONS]")
    print()
    print("Options:")
    print("  -v, --verbose           Show detailed output for each test")
    print(
        "  --planets-only          Test only planetary positions (not ayanamsha values)"
    )
    print("  --ayanamsha-only        Test only ayanamsha values (not positions)")
    print("  -h, --help              Show this help message")
    print()


def main():
    """Main entry point."""
    args = parse_args(sys.argv)

    if args["help"]:
        print_help()
        sys.exit(0)

    planets_only = "--planets-only" in sys.argv
    ayanamsha_only = "--ayanamsha-only" in sys.argv

    passed, total = run_all_comparisons(
        verbose=args["verbose"],
        planets_only=planets_only,
        ayanamsha_only=ayanamsha_only,
    )

    # Exit with appropriate code
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
