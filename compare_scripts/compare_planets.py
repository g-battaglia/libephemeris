"""
Comprehensive Planetary Calculations Comparison Script.

Validates all planetary calculations between SwissEphemeris and libephemeris:
- All major planets (Sun through Pluto)
- Multiple calculation modes (geocentric, heliocentric, barycentric)
- Position (longitude, latitude, distance) and velocity
- Various time periods and edge cases
"""

import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
import sys

sys.path.insert(0, "/Users/giacomo/dev/libephemeris/compare_scripts")
from comparison_utils import *


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets to test
PLANETS = {
    SE_SUN: "Sun",
    SE_MOON: "Moon",
    SE_MERCURY: "Mercury",
    SE_VENUS: "Venus",
    SE_MARS: "Mars",
    SE_JUPITER: "Jupiter",
    SE_SATURN: "Saturn",
    SE_URANUS: "Uranus",
    SE_NEPTUNE: "Neptune",
    SE_PLUTO: "Pluto",
}

# Calculation flags to test
CALC_MODES = {
    SEFLG_SWIEPH: "Geocentric",
    SEFLG_SWIEPH | SEFLG_SPEED: "Geocentric+Speed",
    SEFLG_SWIEPH | SEFLG_HELCTR: "Heliocentric",
    SEFLG_SWIEPH | SEFLG_BARYCTR: "Barycentric",
}

# Test dates (various epochs and edge cases)
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1900, 1, 1, 0.0, "Early 1900s"),
    (2100, 12, 31, 23.999, "Late 2100s"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_planet_position(
    planet: int,
    planet_name: str,
    jd: float,
    flag: int,
    mode_name: str,
    date_str: str,
    verbose: bool = False,
) -> tuple[bool, float, bool]:
    """
    Compare planetary position calculation.
    
    Returns:
        (passed, max_diff, error_occurred)
    """
    try:
        # SwissEphemeris
        pos_swe, ret_swe = swe.calc_ut(jd, planet, flag)
        
        # LibEphemeris
        pos_py, ret_py = pyephem.swe_calc_ut(jd, planet, flag)
        
    except Exception as e:
        if verbose:
            print(f"[{date_str}] [{planet_name}] [{mode_name}]: ERROR {e}")
        return False, 0.0, True
    
    # Compare coordinates
    diff_lon = angular_diff(pos_swe[0], pos_py[0])
    diff_lat = abs(pos_swe[1] - pos_py[1])
    diff_dist = abs(pos_swe[2] - pos_py[2])
    
    # Determine tolerance based on position
    # Geocentric positions should be very precise
    # Heliocentric/barycentric can have slightly larger differences
    if flag & SEFLG_HELCTR or flag & SEFLG_BARYCTR:
        lon_tol = Tolerances.LONGITUDE_RELAXED
        lat_tol = Tolerances.LATITUDE_RELAXED
        dist_tol = Tolerances.DISTANCE_RELAXED
    else:
        lon_tol = Tolerances.LONGITUDE_STRICT
        lat_tol = Tolerances.LATITUDE_STRICT
        dist_tol = Tolerances.DISTANCE_STRICT
    
    # Check if within tolerances
    passed = (
        diff_lon < lon_tol
        and diff_lat < lat_tol
        and diff_dist < dist_tol
    )
    
    max_diff = max(diff_lon, diff_lat * 100)  # Scale lat for comparison
    
    # Compare velocities if SPEED flag is set
    if flag & SEFLG_SPEED:
        diff_lon_speed = abs(pos_swe[3] - pos_py[3])
        diff_lat_speed = abs(pos_swe[4] - pos_py[4])
        diff_dist_speed = abs(pos_swe[5] - pos_py[5])
        
        vel_passed = (
            diff_lon_speed < Tolerances.VELOCITY_ANGULAR
            and diff_lat_speed < Tolerances.VELOCITY_ANGULAR
            and diff_dist_speed < Tolerances.VELOCITY_RADIAL
        )
        
        passed = passed and vel_passed
        max_diff = max(max_diff, diff_lon_speed * 100)
    
    if verbose:
        print(f"\n{'=' * 80}")
        print(f"{planet_name} - {date_str} - {mode_name}")
        print(f"{'=' * 80}")
        print("\nPosition:")
        print(f"  Longitude:  SWE={pos_swe[0]:.8f}°  PY={pos_py[0]:.8f}°  "
              f"Diff={diff_lon:.8f}° {format_status(diff_lon < lon_tol)}")
        print(f"  Latitude:   SWE={pos_swe[1]:.8f}°  PY={pos_py[1]:.8f}°  "
              f"Diff={diff_lat:.8f}° {format_status(diff_lat < lat_tol)}")
        print(f"  Distance:   SWE={pos_swe[2]:.8f} AU  PY={pos_py[2]:.8f} AU  "
              f"Diff={diff_dist:.8f} {format_status(diff_dist < dist_tol)}")
        
        if flag & SEFLG_SPEED:
            print("\nVelocity:")
            print(f"  dLon/dt:    SWE={pos_swe[3]:.8f}°/d  PY={pos_py[3]:.8f}°/d  "
                  f"Diff={diff_lon_speed:.8f}")
            print(f"  dLat/dt:    SWE={pos_swe[4]:.8f}°/d  PY={pos_py[4]:.8f}°/d  "
                  f"Diff={diff_lat_speed:.8f}")
            print(f"  dDist/dt:   SWE={pos_swe[5]:.8f} AU/d  PY={pos_py[5]:.8f} AU/d  "
                  f"Diff={diff_dist_speed:.8f}")
        
        print(f"\nStatus: {'PASSED ✓' if passed else 'FAILED ✗'}")
    else:
        # Compact one-line format
        status = format_status(passed)
        print(f"[{date_str}] [{planet_name:10}] [{mode_name:20}] "
              f"Lon={pos_swe[0]:9.4f}/{pos_py[0]:9.4f} "
              f"Lat={pos_swe[1]:7.4f}/{pos_py[1]:7.4f} "
              f"Dist={pos_swe[2]:8.5f}/{pos_py[2]:8.5f} "
              f"MaxDiff={format_diff(max_diff, 6, 8)} {status}")
    
    return passed, max_diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose: bool = False, planet_filter: str = "all") -> tuple:
    """
    Run all planetary comparison tests.
    
    Args:
        verbose: If True, print detailed output
        planet_filter: 'all', 'inner' (Sun-Mars), or 'outer' (Jupiter-Pluto)
        
    Returns:
        (passed_count, total_count)
    """
    print_header("PLANETARY CALCULATIONS COMPARISON")
    
    # Filter planets
    if planet_filter == "inner":
        planets = {k: v for k, v in PLANETS.items() if k <= SE_MARS}
    elif planet_filter == "outer":
        planets = {k: v for k, v in PLANETS.items() if k >= SE_JUPITER}
    else:
        planets = PLANETS
    
    stats = TestStatistics()
    
    for year, month, day, hour, date_str in TEST_DATES:
        jd = swe.julday(year, month, day, hour)
        
        if not verbose:
            print(f"\n{'=' * 80}")
            print(f"DATE: {date_str} ({year}-{month:02d}-{day:02d})")
            print(f"{'=' * 80}")
        
        for planet, planet_name in planets.items():
            for flag, mode_name in CALC_MODES.items():
                passed, diff, error = compare_planet_position(
                    planet=planet,
                    planet_name=planet_name,
                    jd=jd,
                    flag=flag,
                    mode_name=mode_name,
                    date_str=date_str,
                    verbose=verbose,
                )
                
                stats.add_result(passed, diff, error)
    
    # Print summary
    stats.print_summary("PLANETARY CALCULATIONS SUMMARY")
    
    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_planets.py [OPTIONS]")
    print()
    print("Options:")
    print("  -v, --verbose       Show detailed output for each test")
    print("  -q, --quiet         Suppress all output except summary")
    print("  --inner             Test only inner planets (Sun-Mars)")
    print("  --outer             Test only outer planets (Jupiter-Pluto)")
    print("  -h, --help          Show this help message")
    print()


def main():
    """Main entry point."""
    args = parse_args(sys.argv)
    
    if args["help"]:
        print_help()
        sys.exit(0)
    
    # Determine planet filter
    if "--inner" in sys.argv:
        planet_filter = "inner"
    elif "--outer" in sys.argv:
        planet_filter = "outer"
    else:
        planet_filter = "all"
    
    passed, total = run_all_comparisons(
        verbose=args["verbose"],
        planet_filter=planet_filter,
    )
    
    # Exit with appropriate code
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
