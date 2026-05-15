#!/usr/bin/env python3
import os

"""
Validation script for TNO position precision after Uranus/Neptune perturbations.

Compares Eris, Makemake, Ixion, and Orcus positions against pyswisseph over
a 50-year span (2000-2050) to document precision improvements from adding
Uranus and Neptune perturbations to the secular theory.

This script validates that the Keplerian + secular perturbation model in
libephemeris provides reasonable accuracy for astrological applications,
while documenting the expected precision limitations compared to full
numerical integration (as used in Swiss Ephemeris).

Usage:
    python validate_tno_precision.py [--verbose] [--json]
"""

import sys
import json
import math
from dataclasses import dataclass, asdict
from typing import Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import ERIS, MAKEMAKE, IXION, ORCUS
from libephemeris.minor_bodies import MINOR_BODY_ELEMENTS
from comparison_utils import angular_diff, format_status


# TNO definitions: (libephemeris_id, swisseph_ast_number, name, resonance_info)
TNOS = [
    (ERIS, 136199, "Eris", "scattered disk, ~68 AU"),
    (MAKEMAKE, 136472, "Makemake", "classical KBO, ~45 AU"),
    (IXION, 28978, "Ixion", "plutino (2:3 resonance), ~39 AU"),
    (ORCUS, 90482, "Orcus", "plutino (2:3 resonance), ~39 AU"),
]

# Time span: 2000-2050 (50 years)
# Using annual samples to cover the full span
START_YEAR = 2000
END_YEAR = 2050
SAMPLE_INTERVAL_YEARS = 1  # Sample every year


@dataclass
class TNOValidationResult:
    """Results for a single TNO validation over the time span."""

    name: str
    resonance_info: str
    body_id: int

    # Statistics over the time span
    num_samples: int = 0
    num_skipped: int = 0  # Due to SwissEph data file issues

    # Longitude differences (degrees)
    lon_diff_max: float = 0.0
    lon_diff_min: float = float("inf")
    lon_diff_mean: float = 0.0
    lon_diff_rms: float = 0.0

    # Latitude differences (degrees)
    lat_diff_max: float = 0.0
    lat_diff_min: float = float("inf")
    lat_diff_mean: float = 0.0

    # Distance differences (AU)
    dist_diff_max: float = 0.0
    dist_diff_min: float = float("inf")
    dist_diff_mean: float = 0.0

    # Summary
    passed: bool = False
    tolerance_lon: float = 10.0  # degrees - relaxed for Keplerian approximation
    tolerance_lat: float = 5.0  # degrees

    # Worst case info
    worst_jd: float = 0.0
    worst_year: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


@dataclass
class ValidationSummary:
    """Overall validation summary for all TNOs."""

    total_bodies: int
    total_samples: int
    total_skipped: int
    all_passed: bool

    # Aggregate statistics
    overall_lon_max: float
    overall_lon_mean: float
    overall_lat_max: float
    overall_dist_max: float

    results: list  # List of TNOValidationResult dicts

    # Precision assessment
    precision_grade: str  # "excellent", "good", "acceptable", "poor"
    precision_note: str


def calc_swe_tno_position(jd: float, ast_num: int) -> Optional[tuple]:
    """
    Calculate TNO position using Swiss Ephemeris.

    Args:
        jd: Julian Day (UT)
        ast_num: Asteroid catalog number (e.g., 136199 for Eris)

    Returns:
        (lon, lat, dist) tuple or None if data file missing
    """
    # Swiss Ephemeris uses AST_OFFSET + asteroid_number for asteroids
    swe_body = swe.AST_OFFSET + ast_num

    try:
        pos, _ = swe.calc_ut(jd, swe_body, 0)
        return pos[0], pos[1], pos[2]
    except swe.Error:
        # Data file not available
        return None


def calc_libeph_tno_position(jd: float, body_id: int) -> tuple:
    """
    Calculate TNO position using libephemeris.

    Args:
        jd: Julian Day (UT)
        body_id: libephemeris body ID (e.g., ERIS)

    Returns:
        (lon, lat, dist) tuple
    """
    pos, _ = ephem.calc_ut(jd, body_id, 0)
    return pos[0], pos[1], pos[2]


def validate_tno(
    body_id: int, ast_num: int, name: str, resonance_info: str, verbose: bool = False
) -> TNOValidationResult:
    """
    Validate a single TNO against Swiss Ephemeris over 2000-2050.

    Args:
        body_id: libephemeris body ID
        ast_num: Swiss Ephemeris asteroid number
        name: Body name
        resonance_info: Description of orbital resonance
        verbose: Print detailed output

    Returns:
        TNOValidationResult with statistics
    """
    result = TNOValidationResult(
        name=name, resonance_info=resonance_info, body_id=body_id
    )

    lon_diffs = []
    lat_diffs = []
    dist_diffs = []

    if verbose:
        print(f"\n{'=' * 60}")
        print(f"Validating {name} ({resonance_info})")
        print(f"{'=' * 60}")

    for year in range(START_YEAR, END_YEAR + 1, SAMPLE_INTERVAL_YEARS):
        # Use January 1 at noon for each year
        jd = swe.julday(year, 1, 1, 12.0)

        # Get Swiss Ephemeris position (reference)
        swe_pos = calc_swe_tno_position(jd, ast_num)

        if swe_pos is None:
            result.num_skipped += 1
            if verbose and result.num_skipped == 1:
                print(f"  WARNING: SwissEph data file not available for {name}")
            continue

        # Get libephemeris position
        lib_pos = calc_libeph_tno_position(jd, body_id)

        # Calculate differences
        diff_lon = angular_diff(swe_pos[0], lib_pos[0])
        diff_lat = abs(swe_pos[1] - lib_pos[1])
        diff_dist = abs(swe_pos[2] - lib_pos[2])

        lon_diffs.append(diff_lon)
        lat_diffs.append(diff_lat)
        dist_diffs.append(diff_dist)

        result.num_samples += 1

        # Track worst case
        if diff_lon > result.lon_diff_max:
            result.lon_diff_max = diff_lon
            result.worst_jd = jd
            result.worst_year = year

        result.lon_diff_min = min(result.lon_diff_min, diff_lon)
        result.lat_diff_max = max(result.lat_diff_max, diff_lat)
        result.lat_diff_min = min(result.lat_diff_min, diff_lat)
        result.dist_diff_max = max(result.dist_diff_max, diff_dist)
        result.dist_diff_min = min(result.dist_diff_min, diff_dist)

        if verbose and year % 10 == 0:
            status = "OK" if diff_lon < result.tolerance_lon else "HIGH"
            print(
                f"  {year}: lon_diff={diff_lon:7.3f}° lat_diff={diff_lat:6.3f}° "
                f"dist_diff={diff_dist:.4f} AU [{status}]"
            )

    # Calculate statistics
    if result.num_samples > 0:
        result.lon_diff_mean = sum(lon_diffs) / len(lon_diffs)
        result.lat_diff_mean = sum(lat_diffs) / len(lat_diffs)
        result.dist_diff_mean = sum(dist_diffs) / len(dist_diffs)

        # RMS for longitude
        result.lon_diff_rms = math.sqrt(sum(d * d for d in lon_diffs) / len(lon_diffs))

        # Check pass/fail
        result.passed = (
            result.lon_diff_max < result.tolerance_lon
            and result.lat_diff_max < result.tolerance_lat
        )
    else:
        result.passed = False  # No data to validate
        result.lon_diff_min = 0.0
        result.lat_diff_min = 0.0
        result.dist_diff_min = 0.0

    if verbose:
        print(f"\n  Summary for {name}:")
        print(f"    Samples: {result.num_samples}, Skipped: {result.num_skipped}")
        print(
            f"    Longitude: max={result.lon_diff_max:.3f}°, "
            f"mean={result.lon_diff_mean:.3f}°, RMS={result.lon_diff_rms:.3f}°"
        )
        print(
            f"    Latitude:  max={result.lat_diff_max:.3f}°, "
            f"mean={result.lat_diff_mean:.3f}°"
        )
        print(
            f"    Distance:  max={result.dist_diff_max:.4f} AU, "
            f"mean={result.dist_diff_mean:.4f} AU"
        )
        print(f"    Worst case: year {result.worst_year}")
        print(f"    Status: {format_status(result.passed)}")

    return result


def validate_all_tnos(verbose: bool = False) -> ValidationSummary:
    """
    Validate all TNOs against Swiss Ephemeris.

    Args:
        verbose: Print detailed output

    Returns:
        ValidationSummary with aggregate statistics
    """
    results = []
    total_samples = 0
    total_skipped = 0

    for body_id, ast_num, name, resonance_info in TNOS:
        result = validate_tno(body_id, ast_num, name, resonance_info, verbose)
        results.append(result)
        total_samples += result.num_samples
        total_skipped += result.num_skipped

    # Aggregate statistics (only from results with valid samples)
    valid_results = [r for r in results if r.num_samples > 0]

    if valid_results:
        overall_lon_max = max(r.lon_diff_max for r in valid_results)
        overall_lon_mean = sum(r.lon_diff_mean for r in valid_results) / len(
            valid_results
        )
        overall_lat_max = max(r.lat_diff_max for r in valid_results)
        overall_dist_max = max(r.dist_diff_max for r in valid_results)
        all_passed = all(r.passed for r in valid_results)
    else:
        overall_lon_max = 0.0
        overall_lon_mean = 0.0
        overall_lat_max = 0.0
        overall_dist_max = 0.0
        all_passed = False

    # Assess precision grade
    if overall_lon_max < 1.0:
        precision_grade = "excellent"
        precision_note = (
            "Sub-degree precision achieved. Suitable for most astrological "
            "applications including aspect calculations."
        )
    elif overall_lon_max < 3.0:
        precision_grade = "good"
        precision_note = (
            "Good precision (< 3°). Suitable for general astrological use, "
            "sign positions may need verification near boundaries."
        )
    elif overall_lon_max < 10.0:
        precision_grade = "acceptable"
        precision_note = (
            "Acceptable precision for astrological applications. "
            "Sign positions are reliable, close aspects should be verified."
        )
    else:
        precision_grade = "poor"
        precision_note = (
            "Precision exceeds 10°. Results should be treated as approximate. "
            "Consider using SPK kernels for better accuracy."
        )

    summary = ValidationSummary(
        total_bodies=len(TNOS),
        total_samples=total_samples,
        total_skipped=total_skipped,
        all_passed=all_passed,
        overall_lon_max=overall_lon_max,
        overall_lon_mean=overall_lon_mean,
        overall_lat_max=overall_lat_max,
        overall_dist_max=overall_dist_max,
        results=[r.to_dict() for r in results],
        precision_grade=precision_grade,
        precision_note=precision_note,
    )

    return summary


def print_summary(summary: ValidationSummary):
    """Print formatted validation summary."""
    print("\n" + "=" * 80)
    print("TNO PRECISION VALIDATION SUMMARY")
    print("Comparing libephemeris (Keplerian + secular perturbations)")
    print("against Swiss Ephemeris (full numerical integration)")
    print("Time span: 2000-2050 (50 years)")
    print("=" * 80)

    print(f"\nBodies validated: {summary.total_bodies}")
    print(f"Total samples: {summary.total_samples}")
    print(f"Skipped (no SWE data): {summary.total_skipped}")
    print(f"Overall status: {format_status(summary.all_passed)}")

    print("\n" + "-" * 60)
    print("PRECISION STATISTICS")
    print("-" * 60)
    print(f"  Maximum longitude error:  {summary.overall_lon_max:.3f}°")
    print(f"  Mean longitude error:     {summary.overall_lon_mean:.3f}°")
    print(f"  Maximum latitude error:   {summary.overall_lat_max:.3f}°")
    print(f"  Maximum distance error:   {summary.overall_dist_max:.4f} AU")

    print("\n" + "-" * 60)
    print("PRECISION ASSESSMENT")
    print("-" * 60)
    print(f"  Grade: {summary.precision_grade.upper()}")
    print(f"  {summary.precision_note}")

    print("\n" + "-" * 60)
    print("PER-BODY RESULTS")
    print("-" * 60)

    for r in summary.results:
        if r["num_samples"] > 0:
            print(
                f"  {r['name']:12s}: max_lon={r['lon_diff_max']:6.2f}° "
                f"mean_lon={r['lon_diff_mean']:5.2f}° "
                f"max_lat={r['lat_diff_max']:5.2f}° "
                f"{format_status(r['passed'])}"
            )
        else:
            print(f"  {r['name']:12s}: SKIPPED (no SwissEph data)")

    print("\n" + "-" * 60)
    print("NOTES")
    print("-" * 60)
    print("  - Tolerances: longitude < 10°, latitude < 5°")
    print("  - Keplerian + secular perturbation model has inherent limitations")
    print("  - Plutinos (Ixion, Orcus) are in Neptune resonance - special case")
    print("  - For research-grade precision, use SPK kernels or Swiss Ephemeris")
    print("=" * 80)


def main():
    """Main entry point."""
    args = sys.argv[1:]
    verbose = "--verbose" in args or "-v" in args
    output_json = "--json" in args

    if "--help" in args or "-h" in args:
        print(__doc__)
        return 0

    print("Validating TNO positions against Swiss Ephemeris...")
    print("(Eris, Makemake, Ixion, Orcus over 2000-2050)")

    summary = validate_all_tnos(verbose=verbose)

    if output_json:
        print(json.dumps(asdict(summary), indent=2))
    else:
        print_summary(summary)

    return 0 if summary.all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
