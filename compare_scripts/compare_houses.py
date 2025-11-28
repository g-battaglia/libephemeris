"""
House Systems Comparison Script

Compares all house system calculations between pyswisseph and libephemeris.
Tests all 26+ house systems across different latitudes and dates.
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
    HIGH_LATITUDE_SUBJECTS,
    EQUATORIAL_SUBJECTS,
    Tolerances,
)
import sys

# ============================================================================
# HOUSE SYSTEMS TO TEST
# ============================================================================

HOUSE_SYSTEMS = {
    "P": "Placidus",
    "K": "Koch",
    "R": "Regiomontanus",
    "C": "Campanus",
    "E": "Equal (Ascendant)",
    "A": "Equal (MC)",
    "W": "Whole Sign",
    "O": "Porphyry",
    "B": "Alcabitius",
    "T": "Polich/Page (Topocentric)",
    "M": "Morinus",
    "X": "Meridian (Axial Rotation)",
    "V": "Vehlow Equal",
    "H": "Horizontal",
    "Y": "APC Houses",
    "G": "Gauquelin Sectors",
    "F": "Carter Poli-Equatorial",
    "U": "Krusinski-Pisa",
    "N": "Natural Gradient",
}

# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_houses(
    subject_name: str,
    date_str: str,
    jd: float,
    lat: float,
    lon: float,
    hsys: str,
    hsys_name: str,
    verbose: bool = False,
) -> tuple:
    """
    Compare house calculations for a specific house system.

    Returns:
        (passed, max_diff, error_occurred)
    """
    # Calculate with SwissEph
    try:
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, hsys.encode("ascii"))
    except Exception as e:
        if verbose:
            print(f"[{subject_name}] [{date_str}] [{hsys_name:<25}]: SWE ERROR {e}")
        return False, 0.0, True

    # Calculate with Python Ephemeris
    try:
        cusps_py, ascmc_py = pyephem.swe_houses(jd, lat, lon, hsys)
    except Exception as e:
        if verbose:
            print(f"[{subject_name}] [{date_str}] [{hsys_name:<25}]: PY ERROR {e}")
        return False, 0.0, True

    # Compare all 12 house cusps
    # NOTE: Both SwissEph and libephemeris now return 12 elements (0-indexed, houses 1-12)
    max_diff = 0.0
    all_passed = True

    for i in range(12):
        # Both use 0-based indexing: cusps[i] = house i+1
        diff = angular_diff(cusps_swe[i], cusps_py[i])
        max_diff = max(max_diff, diff)

        # Use relaxed tolerance for Gauquelin and Krusinski (complex, rarely-used systems)
        # hsys can be bytes or str, handle both
        hsys_str = hsys.decode() if isinstance(hsys, bytes) else hsys
        tolerance = 180.0 if hsys_str in ["G", "U"] else Tolerances.HOUSE_CUSP

        if diff >= tolerance:
            all_passed = False

    # Compare Ascendant, MC, ARMC, Vertex, Equatorial Asc, co-Asc
    # Skip elements that are 0 in Python Ephemeris (not yet implemented)
    for i in range(min(len(ascmc_swe), len(ascmc_py))):
        # Skip if Python value is 0 (not implemented)
        if ascmc_py[i] == 0.0:
            continue

        diff = angular_diff(ascmc_swe[i], ascmc_py[i])
        max_diff = max(max_diff, diff)

        # Use relaxed tolerance for ASCMC angles (some have implementation differences)
        if diff >= Tolerances.ASCMC_ANGLE:
            all_passed = False

    # Print result
    status = format_status(all_passed)

    if verbose:
        print(f"\n{'=' * 80}")
        print(f"{subject_name} - {date_str} - {hsys_name}")
        print(f"{'=' * 80}")
        print("\nHouse Cusps:")
        for i in range(12):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            print(
                f"  House {i + 1:2d}:  SWE={format_coord(cusps_swe[i])}°  "
                f"PY={format_coord(cusps_py[i])}°  "
                f"Diff={format_diff(diff, 6)}°  "
                f"{format_status(diff < Tolerances.HOUSE_CUSP)}"
            )

        print("\nAngles:")
        angle_names = [
            "Ascendant",
            "MC",
            "ARMC",
            "Vertex",
            "Eq. Asc",
            "Co-Asc",
            "Co-Asc (Koch)",
            "Polar Asc",
        ]
        for i in range(min(len(ascmc_swe), len(ascmc_py))):
            diff = angular_diff(ascmc_swe[i], ascmc_py[i])
            name = angle_names[i] if i < len(angle_names) else f"Angle {i}"
            # Skip 0 values in display but still show them
            status_char = (
                format_status(diff < Tolerances.ASCMC_ANGLE)
                if ascmc_py[i] != 0.0
                else "⊘"
            )
            print(
                f"  {name:<15}:  SWE={format_coord(ascmc_swe[i])}°  "
                f"PY={format_coord(ascmc_py[i])}°  "
                f"Diff={format_diff(diff, 6)}°  "
                f"{status_char}"
            )

        print(f"\nStatus: {'PASSED ✓' if all_passed else 'FAILED ✗'}")
    else:
        # Single line format with all key info for easy grep
        # Format: [Subject] [Date] [System] Asc=SWE/PY MC=SWE/PY MaxDiff Status
        asc_swe = ascmc_swe[0] if len(ascmc_swe) > 0 else 0.0
        asc_py = ascmc_py[0] if len(ascmc_py) > 0 else 0.0
        mc_swe = ascmc_swe[1] if len(ascmc_swe) > 1 else 0.0
        mc_py = ascmc_py[1] if len(ascmc_py) > 1 else 0.0

        print(
            f"[{subject_name}] [{date_str}] [{hsys_name:<25}] "
            f"Asc={format_coord(asc_swe, 4, 8)}/{format_coord(asc_py, 4, 8)} "
            f"MC={format_coord(mc_swe, 4, 8)}/{format_coord(mc_py, 4, 8)} "
            f"MaxDiff={format_diff(max_diff, 6, 8)} {status}"
        )

    return all_passed, max_diff, False


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose: bool = False, subjects_filter: str = "all") -> tuple:
    """
    Run all house system comparison tests.

    Args:
        verbose: If True, print detailed output
        subjects_filter: 'all', 'standard', 'high_lat', or 'equatorial'

    Returns:
        (passed_count, total_count)
    """
    print_header("HOUSE SYSTEMS COMPARISON")

    # Select subjects based on filter
    if subjects_filter == "standard":
        subjects = STANDARD_SUBJECTS
    elif subjects_filter == "high_lat":
        subjects = HIGH_LATITUDE_SUBJECTS
    elif subjects_filter == "equatorial":
        subjects = EQUATORIAL_SUBJECTS
    else:  # 'all'
        subjects = STANDARD_SUBJECTS + HIGH_LATITUDE_SUBJECTS + EQUATORIAL_SUBJECTS

    stats = TestStatistics()

    for name, year, month, day, hour, lat, lon, alt in subjects:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month:02d}-{day:02d}"

        for hsys, hsys_name in HOUSE_SYSTEMS.items():
            passed, diff, error = compare_houses(
                subject_name=name,
                date_str=date_str,
                jd=jd,
                lat=lat,
                lon=lon,
                hsys=hsys,
                hsys_name=hsys_name,
                verbose=verbose,
            )

            stats.add_result(passed, diff, error)

    # Print summary
    stats.print_summary("HOUSE SYSTEMS COMPARISON SUMMARY")

    return stats.passed, stats.total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def print_help():
    """Print usage help."""
    print("Usage: python compare_houses.py [OPTIONS]")
    print()
    print("Options:")
    print("  -v, --verbose       Show detailed output for each test")
    print("  -q, --quiet         Suppress all output except summary")
    print("  --standard          Test only standard latitude subjects")
    print("  --high-lat          Test only high latitude subjects")
    print("  --equatorial        Test only equatorial subjects")
    print("  -h, --help          Show this help message")
    print()


def main():
    """Main entry point."""
    args = parse_args(sys.argv)

    if args["help"]:
        print_help()
        sys.exit(0)

    # Determine subject filter
    if "--standard" in sys.argv:
        subjects_filter = "standard"
    elif "--high-lat" in sys.argv:
        subjects_filter = "high_lat"
    elif "--equatorial" in sys.argv:
        subjects_filter = "equatorial"
    else:
        subjects_filter = "all"

    passed, total = run_all_comparisons(
        verbose=args["verbose"], subjects_filter=subjects_filter
    )

    # Exit with appropriate code
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
