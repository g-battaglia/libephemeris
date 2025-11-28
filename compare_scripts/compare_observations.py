"""
Comprehensive comparison module for libephemeris vs pyswisseph.

This module provides detailed comparisons across all observation modes,
coordinate systems, and calculated values (position, velocity, declination).
"""

import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import *
from typing import Tuple, List, Optional
import sys

# ============================================================================
# TEST SUBJECTS
# ============================================================================

SUBJECTS = [
    # (Name, Year, Month, Day, Hour, Lat, Lon, Alt)
    ("Standard J2000", 2000, 1, 1, 12.0, 0.0, 0.0, 0),
    ("Einstein (Ulm)", 1879, 3, 14, 11.5, 48.4011, 9.9876, 478),
    ("Gandhi (Porbandar)", 1869, 10, 2, 7.2, 21.6417, 69.6293, 0),
    ("Mandela (Mvezo)", 1918, 7, 18, 14.0, -31.9566, 28.5133, 0),
    ("Tromso (High Lat)", 1990, 1, 15, 12.0, 69.6492, 18.9553, 0),
    ("McMurdo (High Lat)", 2005, 6, 21, 0.0, -77.8463, 166.6681, 0),
]

# ============================================================================
# COMPARISON CONFIGURATIONS
# ============================================================================

# Format: (mode_name, planet_id, planet_name, flags)
COMPARISON_MODES = [
    ("Topocentric", SE_MOON, "Moon", SEFLG_TOPOCTR | SEFLG_SPEED),
    ("Heliocentric", SE_MARS, "Mars", SEFLG_HELCTR | SEFLG_SPEED),
    ("Barycentric", SE_JUPITER, "Jupiter", SEFLG_BARYCTR | SEFLG_SPEED),
    ("TruePos", SE_VENUS, "Venus", SEFLG_TRUEPOS | SEFLG_SPEED),
    ("J2000 Equ", SE_SUN, "Sun", SEFLG_J2000 | SEFLG_EQUATORIAL | SEFLG_SPEED),
    ("J2000 Ecl", SE_SATURN, "Saturn", SEFLG_J2000 | SEFLG_SPEED),
]

# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class Tolerances:
    """Tolerance thresholds for different comparison types."""

    # Position tolerances (degrees)
    LONGITUDE_STRICT = 0.001  # For geocentric/topocentric
    LONGITUDE_RELAXED = 0.03  # For heliocentric/barycentric
    LATITUDE_STRICT = 0.001  # For declination/latitude
    LATITUDE_RELAXED = 0.03

    # Distance tolerances (AU)
    DISTANCE_STRICT = 0.0001
    DISTANCE_RELAXED = 0.01

    # Velocity tolerances (degrees/day or AU/day)
    VELOCITY_ANGULAR = 0.01  # degrees/day
    VELOCITY_RADIAL = 0.001  # AU/day


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360° wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def format_coord(value: float, decimals: int = 6) -> str:
    """Format coordinate value with consistent width."""
    return f"{value:10.{decimals}f}"


def format_diff(value: float, decimals: int = 8) -> str:
    """Format difference value with consistent width."""
    return f"{value:10.{decimals}f}"


# ============================================================================
# COMPARISON RESULT CLASS
# ============================================================================


class ComparisonResult:
    """Stores and formats comparison results."""

    def __init__(
        self, subject_name: str, date_str: str, mode_name: str, planet_name: str
    ):
        self.subject_name = subject_name
        self.date_str = date_str
        self.mode_name = mode_name
        self.planet_name = planet_name

        # SwissEph values
        self.lon_swe = 0.0
        self.lat_swe = 0.0
        self.dist_swe = 0.0
        self.lon_speed_swe = 0.0
        self.lat_speed_swe = 0.0
        self.dist_speed_swe = 0.0

        # Python Ephemeris values
        self.lon_py = 0.0
        self.lat_py = 0.0
        self.dist_py = 0.0
        self.lon_speed_py = 0.0
        self.lat_speed_py = 0.0
        self.dist_speed_py = 0.0

        # Status
        self.error_swe: Optional[str] = None
        self.error_py: Optional[str] = None
        self.passed = False

    def calculate_diffs(self):
        """Calculate all differences."""
        self.diff_lon = angular_diff(self.lon_swe, self.lon_py)
        self.diff_lat = abs(self.lat_swe - self.lat_py)
        self.diff_dist = abs(self.dist_swe - self.dist_py)
        self.diff_lon_speed = abs(self.lon_speed_swe - self.lon_speed_py)
        self.diff_lat_speed = abs(self.lat_speed_swe - self.lat_speed_py)
        self.diff_dist_speed = abs(self.dist_speed_swe - self.dist_speed_py)

    def check_passed(self, use_relaxed: bool = False) -> bool:
        """Check if comparison passed within tolerances."""
        if self.error_swe or self.error_py:
            self.passed = False
            return False

        lon_tol = (
            Tolerances.LONGITUDE_RELAXED if use_relaxed else Tolerances.LONGITUDE_STRICT
        )
        lat_tol = (
            Tolerances.LATITUDE_RELAXED if use_relaxed else Tolerances.LATITUDE_STRICT
        )
        dist_tol = (
            Tolerances.DISTANCE_RELAXED if use_relaxed else Tolerances.DISTANCE_STRICT
        )

        self.passed = (
            self.diff_lon < lon_tol
            and self.diff_lat < lat_tol
            and self.diff_dist < dist_tol
            and self.diff_lon_speed < Tolerances.VELOCITY_ANGULAR
            and self.diff_lat_speed < Tolerances.VELOCITY_ANGULAR
            and self.diff_dist_speed < Tolerances.VELOCITY_RADIAL
        )
        return self.passed

    def format_summary(self) -> str:
        """Format one-line summary."""
        if self.error_swe:
            return f"[{self.subject_name}] [{self.date_str}] [{self.mode_name:<20}] {self.planet_name:<10}: SWE ERROR {self.error_swe}"
        if self.error_py:
            return f"[{self.subject_name}] [{self.date_str}] [{self.mode_name:<20}] {self.planet_name:<10}: PY ERROR {self.error_py}"

        status = "✓" if self.passed else "✗"
        return (
            f"[{self.subject_name}] [{self.date_str}] [{self.mode_name:<20}] {self.planet_name:<10}: "
            f"Lon {format_diff(self.diff_lon)} | Lat {format_diff(self.diff_lat, 6)} | "
            f"Dist {format_diff(self.diff_dist, 6)} | "
            f"LonSpd {format_diff(self.diff_lon_speed, 6)} | "
            f"LatSpd {format_diff(self.diff_lat_speed, 6)} {status}"
        )

    def format_detailed(self) -> str:
        """Format detailed multi-line output."""
        lines = []
        lines.append(f"\n{'=' * 80}")
        lines.append(
            f"{self.subject_name} - {self.date_str} - {self.mode_name} - {self.planet_name}"
        )
        lines.append(f"{'=' * 80}")

        if self.error_swe:
            lines.append(f"SwissEph ERROR: {self.error_swe}")
            return "\n".join(lines)
        if self.error_py:
            lines.append(f"Python Ephemeris ERROR: {self.error_py}")
            return "\n".join(lines)

        # Position comparison
        lines.append("\nPosition:")
        lines.append(
            f"  Longitude:  SWE={format_coord(self.lon_swe)}°  PY={format_coord(self.lon_py)}°  Diff={format_diff(self.diff_lon)}°"
        )
        lines.append(
            f"  Latitude:   SWE={format_coord(self.lat_swe)}°  PY={format_coord(self.lat_py)}°  Diff={format_diff(self.diff_lat, 6)}°"
        )
        lines.append(
            f"  Distance:   SWE={format_coord(self.dist_swe)} AU  PY={format_coord(self.dist_py)} AU  Diff={format_diff(self.diff_dist, 6)} AU"
        )

        # Velocity comparison
        lines.append("\nVelocity:")
        lines.append(
            f"  Lon Speed:  SWE={format_coord(self.lon_speed_swe, 8)}°/d  PY={format_coord(self.lon_speed_py, 8)}°/d  Diff={format_diff(self.diff_lon_speed, 6)}"
        )
        lines.append(
            f"  Lat Speed:  SWE={format_coord(self.lat_speed_swe, 8)}°/d  PY={format_coord(self.lat_speed_py, 8)}°/d  Diff={format_diff(self.diff_lat_speed, 6)}"
        )
        lines.append(
            f"  Dist Speed: SWE={format_coord(self.dist_speed_swe, 8)} AU/d  PY={format_coord(self.dist_speed_py, 8)} AU/d  Diff={format_diff(self.diff_dist_speed, 6)}"
        )

        status_msg = "PASSED ✓" if self.passed else "FAILED ✗"
        lines.append(f"\nStatus: {status_msg}")

        return "\n".join(lines)


# ============================================================================
# COMPARISON FUNCTIONS
# ============================================================================


def compare_calculation(
    subject_name: str,
    date_str: str,
    jd: float,
    planet_id: int,
    planet_name: str,
    mode_name: str,
    flags: int,
    lat: float = 0,
    lon: float = 0,
    alt: float = 0,
    verbose: bool = False,
) -> ComparisonResult:
    """
    Compare a single planetary calculation between SwissEph and Python Ephemeris.

    Args:
        subject_name: Name of the test subject
        date_str: Date string for display
        jd: Julian day
        planet_id: Planet ID constant
        planet_name: Planet name for display
        mode_name: Mode name for display
        flags: Calculation flags
        lat: Latitude (if topocentric)
        lon: Longitude (if topocentric)
        alt: Altitude (if topocentric)
        verbose: If True, print detailed output

    Returns:
        ComparisonResult object
    """
    result = ComparisonResult(subject_name, date_str, mode_name, planet_name)

    # Set topocentric location if needed
    if flags & SEFLG_TOPOCTR:
        swe.set_topo(lon, lat, alt)
        pyephem.swe_set_topo(lon, lat, alt)

    # Calculate with SwissEph
    try:
        res_swe, _ = swe.calc_ut(jd, planet_id, flags)
        result.lon_swe = res_swe[0]
        result.lat_swe = res_swe[1]
        result.dist_swe = res_swe[2]
        result.lon_speed_swe = res_swe[3]
        result.lat_speed_swe = res_swe[4]
        result.dist_speed_swe = res_swe[5]
    except Exception as e:
        result.error_swe = str(e)
        return result

    # Calculate with Python Ephemeris
    try:
        res_py, _ = pyephem.swe_calc_ut(jd, planet_id, flags)
        result.lon_py = res_py[0]
        result.lat_py = res_py[1]
        result.dist_py = res_py[2]
        result.lon_speed_py = res_py[3]
        result.lat_speed_py = res_py[4]
        result.dist_speed_py = res_py[5]
    except Exception as e:
        result.error_py = str(e)
        return result

    # Calculate differences and check pass/fail
    result.calculate_diffs()

    # Use relaxed tolerances for Heliocentric/Barycentric
    use_relaxed = bool(flags & (SEFLG_HELCTR | SEFLG_BARYCTR))
    result.check_passed(use_relaxed)

    if verbose:
        print(result.format_detailed())
    else:
        print(result.format_summary())

    return result


# ============================================================================
# MAIN COMPARISON RUNNER
# ============================================================================


def run_all_comparisons(verbose: bool = False) -> Tuple[int, int]:
    """
    Run all comparison tests.

    Args:
        verbose: If True, print detailed output for each test

    Returns:
        Tuple of (passed_count, total_count)
    """
    print("=" * 80)
    print("COMPREHENSIVE OBSERVATION MODES COMPARISON")
    print("=" * 80)
    print()

    results: List[ComparisonResult] = []

    for name, year, month, day, hour, lat, lon, alt in SUBJECTS:
        jd = swe.julday(year, month, day, hour)
        date_str = f"{year}-{month}-{day}"

        for mode_name, planet_id, planet_name, flags in COMPARISON_MODES:
            result = compare_calculation(
                subject_name=name,
                date_str=date_str,
                jd=jd,
                planet_id=planet_id,
                planet_name=planet_name,
                mode_name=mode_name,
                flags=flags,
                lat=lat,
                lon=lon,
                alt=alt,
                verbose=verbose,
            )
            results.append(result)

    # Summary statistics
    total = len(results)
    passed = sum(1 for r in results if r.passed)
    failed = sum(
        1 for r in results if not r.passed and not r.error_swe and not r.error_py
    )
    errors = sum(1 for r in results if r.error_swe or r.error_py)

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total tests:   {total}")
    print(f"Passed:        {passed} ✓")
    print(f"Failed:        {failed} ✗")
    print(f"Errors:        {errors}")
    print(f"Pass rate:     {passed / total * 100:.1f}%")
    print("=" * 80)

    return passed, total


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================


def main():
    """Main entry point."""
    verbose = "--verbose" in sys.argv or "-v" in sys.argv
    passed, total = run_all_comparisons(verbose=verbose)

    # Exit with appropriate code
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
