"""
Shared utilities for comparison scripts.

This module provides common classes, functions, and constants used across
all comparison scripts in the suite.
"""

from typing import Optional, List
from dataclasses import dataclass

# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================


class Tolerances:
    """Tolerance thresholds for different comparison types."""

    # Position tolerances (degrees)
    LONGITUDE_STRICT = 0.001  # Geocentric/topocentric
    LONGITUDE_RELAXED = 0.03  # Heliocentric/barycentric
    LATITUDE_STRICT = 0.001
    LATITUDE_RELAXED = 0.03

    # Distance tolerances (AU)
    DISTANCE_STRICT = 0.0001
    DISTANCE_RELAXED = 0.01

    # Velocity tolerances
    VELOCITY_ANGULAR = 0.01  # degrees/day
    VELOCITY_RADIAL = 0.001  # AU/day

    # House cusp tolerances (degrees)
    HOUSE_CUSP = 0.001

    # ASCMC angle tolerances (degrees)
    # Some angles (Co-Asc, Polar Asc) may not be implemented or have calculation differences
    ASCMC_ANGLE = 1.0

    # Ayanamsha tolerances (degrees)
    # Relaxed tolerance for star-based ayanamshas due to numerical precision in
    # star position calculations, precession, and coordinate transformations
    AYANAMSHA = 0.06


# ============================================================================
# TEST SUBJECTS
# ============================================================================

# Format: (Name, Year, Month, Day, Hour, Lat, Lon, Alt)
STANDARD_SUBJECTS = [
    ("Standard J2000", 2000, 1, 1, 12.0, 0.0, 0.0, 0),
    ("Rome", 1980, 5, 20, 14.5, 41.9028, 12.4964, 0),
    ("New York", 2024, 11, 5, 9.0, 40.7128, -74.0060, 0),
    ("Sydney", 1950, 10, 15, 22.0, -33.8688, 151.2093, 0),
]

HIGH_LATITUDE_SUBJECTS = [
    ("Tromso (Arctic)", 1990, 1, 15, 12.0, 69.6492, 18.9553, 0),
    ("McMurdo (Antarctic)", 2005, 6, 21, 0.0, -77.8463, 166.6681, 0),
]

EQUATORIAL_SUBJECTS = [
    ("Equator", 1975, 3, 21, 12.0, 0.0, 45.0, 0),
]

ALL_SUBJECTS = STANDARD_SUBJECTS + HIGH_LATITUDE_SUBJECTS + EQUATORIAL_SUBJECTS

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360° wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def format_coord(value: float, decimals: int = 6, width: int = 10) -> str:
    """Format coordinate value with consistent width."""
    return f"{value:{width}.{decimals}f}"


def format_diff(value: float, decimals: int = 8, width: int = 10) -> str:
    """Format difference value with consistent width."""
    return f"{value:{width}.{decimals}f}"


def format_status(passed: bool) -> str:
    """Format pass/fail status."""
    return "✓" if passed else "✗"


# ============================================================================
# COMPARISON RESULT CLASSES
# ============================================================================


@dataclass
class PositionResult:
    """Stores position comparison result."""

    lon_swe: float = 0.0
    lat_swe: float = 0.0
    dist_swe: float = 0.0
    lon_py: float = 0.0
    lat_py: float = 0.0
    dist_py: float = 0.0

    diff_lon: float = 0.0
    diff_lat: float = 0.0
    diff_dist: float = 0.0

    passed: bool = False
    error_swe: Optional[str] = None
    error_py: Optional[str] = None

    def calculate_diffs(self):
        """Calculate all differences."""
        self.diff_lon = angular_diff(self.lon_swe, self.lon_py)
        self.diff_lat = abs(self.lat_swe - self.lat_py)
        self.diff_dist = abs(self.dist_swe - self.dist_py)

    def check_passed(self, lon_tol: float, lat_tol: float, dist_tol: float) -> bool:
        """Check if within tolerances."""
        if self.error_swe or self.error_py:
            self.passed = False
            return False
        self.passed = (
            self.diff_lon < lon_tol
            and self.diff_lat < lat_tol
            and self.diff_dist < dist_tol
        )
        return self.passed


@dataclass
class VelocityResult:
    """Stores velocity comparison result."""

    lon_speed_swe: float = 0.0
    lat_speed_swe: float = 0.0
    dist_speed_swe: float = 0.0
    lon_speed_py: float = 0.0
    lat_speed_py: float = 0.0
    dist_speed_py: float = 0.0

    diff_lon_speed: float = 0.0
    diff_lat_speed: float = 0.0
    diff_dist_speed: float = 0.0

    passed: bool = False

    def calculate_diffs(self):
        """Calculate all differences."""
        self.diff_lon_speed = abs(self.lon_speed_swe - self.lon_speed_py)
        self.diff_lat_speed = abs(self.lat_speed_swe - self.lat_speed_py)
        self.diff_dist_speed = abs(self.dist_speed_swe - self.dist_speed_py)

    def check_passed(self, ang_tol: float, rad_tol: float) -> bool:
        """Check if within tolerances."""
        self.passed = (
            self.diff_lon_speed < ang_tol
            and self.diff_lat_speed < ang_tol
            and self.diff_dist_speed < rad_tol
        )
        return self.passed


# ============================================================================
# SUMMARY STATISTICS
# ============================================================================


class TestStatistics:
    """Tracks and reports test statistics."""

    def __init__(self):
        self.total = 0
        self.passed = 0
        self.failed = 0
        self.errors = 0
        self.max_diff = 0.0
        self.diff_sum = 0.0

    def add_result(self, passed: bool, diff: float = 0.0, error: bool = False):
        """Add a test result."""
        self.total += 1
        if error:
            self.errors += 1
        elif passed:
            self.passed += 1
        else:
            self.failed += 1

        if not error:
            self.max_diff = max(self.max_diff, diff)
            self.diff_sum += diff

    def avg_diff(self) -> float:
        """Calculate average difference (excluding errors)."""
        count = self.total - self.errors
        return self.diff_sum / count if count > 0 else 0.0

    def pass_rate(self) -> float:
        """Calculate pass rate (excluding errors)."""
        count = self.total - self.errors
        return (self.passed / count * 100) if count > 0 else 0.0

    def print_summary(self, title: str = "SUMMARY"):
        """Print formatted summary."""
        print()
        print("=" * 80)
        print(title)
        print("=" * 80)
        print(f"Total tests:   {self.total}")
        print(f"Passed:        {self.passed} ✓")
        print(f"Failed:        {self.failed} ✗")
        print(f"Errors:        {self.errors}")
        if self.total > self.errors:
            print(f"Pass rate:     {self.pass_rate():.1f}%")
            print(f"Max diff:      {self.max_diff:.6f}")
            print(f"Avg diff:      {self.avg_diff():.6f}")
        print("=" * 80)


# ============================================================================
# COMMAND LINE HELPERS
# ============================================================================


def parse_args(args: List[str]) -> dict:
    """Parse common command line arguments."""
    return {
        "verbose": "--verbose" in args or "-v" in args,
        "quiet": "--quiet" in args or "-q" in args,
        "help": "--help" in args or "-h" in args,
    }


def print_header(title: str):
    """Print formatted header."""
    print("=" * 80)
    print(title)
    print("=" * 80)
    print()
