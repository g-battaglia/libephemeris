"""
Lunar node and Lilith comparison functions.

This module provides comparison utilities for lunar nodes and Lilith
between pyswisseph and libephemeris.
"""

import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    TRUE_NODE,
    MEAN_NODE,
    MEAN_APOG,
    OSCU_APOG,
)


@dataclass
class ComparisonStats:
    """Statistics for comparison runs."""

    total: int = 0
    passed: int = 0
    failed: int = 0
    errors: int = 0
    max_diff: float = 0.0
    diffs: List[float] = field(default_factory=list)

    def pass_rate(self) -> float:
        """Return pass rate as percentage."""
        return (self.passed / self.total * 100) if self.total > 0 else 0.0

    def mean_diff(self) -> float:
        """Return mean difference."""
        return sum(self.diffs) / len(self.diffs) if self.diffs else 0.0


def generate_random_jd(
    start_year: int, end_year: int, count: int, seed: int | None = None
) -> List[Tuple[int, int, int, float, float]]:
    """
    Generate random Julian Day numbers for testing.

    Args:
        start_year: Start year for random dates
        end_year: End year for random dates
        count: Number of dates to generate
        seed: Optional random seed for reproducibility

    Returns:
        List of (year, month, day, hour, jd) tuples
    """
    if seed is not None:
        random.seed(seed)

    dates = []
    for _ in range(count):
        year = random.randint(start_year, end_year)
        month = random.randint(1, 12)
        day = random.randint(1, 28)  # Safe for all months
        hour = random.uniform(0, 24)
        jd = swe.julday(year, month, day, hour)
        dates.append((year, month, day, hour, jd))

    return dates


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def compare_lunar_nodes(
    jd: float, label: str, date_str: str
) -> Dict[str, Tuple[bool, float]]:
    """
    Compare lunar node calculations between implementations.

    Args:
        jd: Julian Day number
        label: Label for the test
        date_str: Date string for logging

    Returns:
        Dictionary with results for mean_node and true_node
    """
    results = {}
    tolerance = 0.3  # degrees - relaxed for true node

    for name, body_id, swe_id in [
        ("mean_node", MEAN_NODE, swe.MEAN_NODE),
        ("true_node", TRUE_NODE, swe.TRUE_NODE),
    ]:
        try:
            pos_swe, _ = swe.calc_ut(jd, swe_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
            diff = angular_diff(pos_swe[0], pos_py[0])
            passed = diff < tolerance
            results[name] = (passed, diff)
        except Exception:
            results[name] = (False, float("inf"))

    return results


def compare_lilith(
    jd: float, label: str, date_str: str
) -> Dict[str, Tuple[bool, float]]:
    """
    Compare Lilith calculations between implementations.

    Args:
        jd: Julian Day number
        label: Label for the test
        date_str: Date string for logging

    Returns:
        Dictionary with results for mean_lilith and true_lilith
    """
    results = {}

    for name, body_id, swe_id, tol in [
        ("mean_lilith", MEAN_APOG, swe.MEAN_APOG, 0.15),
        ("true_lilith", OSCU_APOG, swe.OSCU_APOG, 10.0),
    ]:
        try:
            pos_swe, _ = swe.calc_ut(jd, swe_id, 0)
            pos_py, _ = ephem.calc_ut(jd, body_id, 0)
            diff = angular_diff(pos_swe[0], pos_py[0])
            passed = diff < tol
            results[name] = (passed, diff)
        except Exception:
            results[name] = (False, float("inf"))

    return results


def compare_true_node_precision(
    count: int = 1000, seed: int = 42, tolerance: float = 0.3
) -> ComparisonStats:
    """
    Compare True Node precision across random dates.

    Args:
        count: Number of random dates to test
        seed: Random seed for reproducibility
        tolerance: Maximum acceptable difference in degrees

    Returns:
        ComparisonStats with results
    """
    stats = ComparisonStats()
    dates = generate_random_jd(1950, 2050, count, seed=seed)

    for year, month, day, hour, jd in dates:
        try:
            pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, 0)
            pos_py, _ = ephem.calc_ut(jd, TRUE_NODE, 0)
            diff = angular_diff(pos_swe[0], pos_py[0])

            stats.total += 1
            stats.diffs.append(diff)
            stats.max_diff = max(stats.max_diff, diff)

            if diff < tolerance:
                stats.passed += 1
            else:
                stats.failed += 1
        except Exception:
            stats.errors += 1

    return stats


def compare_mean_node_precision(
    count: int = 1000, seed: int = 42, tolerance: float = 0.01
) -> ComparisonStats:
    """
    Compare Mean Node precision across random dates.

    Args:
        count: Number of random dates to test
        seed: Random seed for reproducibility
        tolerance: Maximum acceptable difference in degrees

    Returns:
        ComparisonStats with results
    """
    stats = ComparisonStats()
    dates = generate_random_jd(1950, 2050, count, seed=seed)

    for year, month, day, hour, jd in dates:
        try:
            pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, 0)
            pos_py, _ = ephem.calc_ut(jd, MEAN_NODE, 0)
            diff = angular_diff(pos_swe[0], pos_py[0])

            stats.total += 1
            stats.diffs.append(diff)
            stats.max_diff = max(stats.max_diff, diff)

            if diff < tolerance:
                stats.passed += 1
            else:
                stats.failed += 1
        except Exception:
            stats.errors += 1

    return stats


def main() -> int:
    """
    Run the comparison and return exit code.

    Returns:
        0 if pass rate >= 80%, 1 otherwise
    """
    print("Comparing True Node precision...")
    stats = compare_true_node_precision()

    print(f"Total tests: {stats.total}")
    print(f"Passed: {stats.passed}")
    print(f"Failed: {stats.failed}")
    print(f"Errors: {stats.errors}")
    print(f"Pass rate: {stats.pass_rate():.1f}%")
    print(f"Max difference: {stats.max_diff:.4f}°")
    print(f"Mean difference: {stats.mean_diff():.4f}°")

    return 0 if stats.pass_rate() >= 80.0 else 1


if __name__ == "__main__":
    exit(main())
