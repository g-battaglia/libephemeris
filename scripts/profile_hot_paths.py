#!/usr/bin/env python3
"""
Performance profiling script for libephemeris.

This script profiles the most commonly used functions to identify hot paths
and optimization opportunities.

Usage:
    python scripts/profile_hot_paths.py

Output:
    - Detailed cProfile statistics
    - Top functions by cumulative time
    - Top functions by total time (self-time)

Provenance:
    Project-authored performance diagnostic using Python's standard ``cProfile``
    and the registered public API. Workload dates and bodies exercise code paths
    but supply no model data. Profiles may justify engineering optimization;
    they are not scientific validation and are never emitted as coefficients.
"""

import cProfile
import pstats
import io
from pstats import SortKey
import time
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    FLG_SPEED,
)


def benchmark_planet_calculations(n_iterations: int = 100) -> None:
    """Benchmark planet position calculations."""
    jd = ephem.julday(2000, 1, 1, 12.0)
    planets = [
        SUN,
        MOON,
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
    ]

    for _ in range(n_iterations):
        for planet in planets:
            ephem.calc_ut(jd, planet, FLG_SPEED)


def benchmark_house_calculations(n_iterations: int = 100) -> None:
    """Benchmark house cusp calculations."""
    jd = ephem.julday(2000, 1, 1, 12.0)
    lat, lon = 41.9, 12.5  # Rome
    house_systems = [
        ord("P"),
        ord("K"),
        ord("R"),
        ord("C"),
        ord("E"),
        ord("W"),
        ord("O"),
    ]

    for _ in range(n_iterations):
        for hsys in house_systems:
            ephem.houses(jd, lat, lon, hsys)


def benchmark_lunar_nodes(n_iterations: int = 100) -> None:
    """Benchmark lunar node calculations."""
    jd = ephem.julday(2000, 1, 1, 12.0)

    for _ in range(n_iterations):
        ephem.calc_ut(jd, MEAN_NODE, FLG_SPEED)
        ephem.calc_ut(jd, TRUE_NODE, FLG_SPEED)
        ephem.calc_ut(jd, MEAN_APOG, FLG_SPEED)


def benchmark_julian_day(n_iterations: int = 1000) -> None:
    """Benchmark Julian Day calculations."""
    for i in range(n_iterations):
        year = 1900 + (i % 200)
        month = (i % 12) + 1
        day = (i % 28) + 1
        hour = (i % 24) + 0.5
        ephem.julday(year, month, day, hour)


def benchmark_ayanamsa(n_iterations: int = 100) -> None:
    """Benchmark ayanamsa calculations."""
    ephem.set_sid_mode(1)  # Lahiri
    jd = ephem.julday(2000, 1, 1, 12.0)

    for _ in range(n_iterations):
        ephem.get_ayanamsa_ut(jd)


def benchmark_full_chart(n_iterations: int = 50) -> None:
    """Benchmark a complete chart calculation (planets + houses + nodes)."""
    jd = ephem.julday(2000, 1, 1, 12.0)
    lat, lon = 41.9, 12.5  # Rome
    planets = [
        SUN,
        MOON,
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
        MEAN_NODE,
        TRUE_NODE,
        MEAN_APOG,
    ]

    for _ in range(n_iterations):
        # Calculate all planets
        for planet in planets:
            ephem.calc_ut(jd, planet, FLG_SPEED)

        # Calculate houses (Placidus)
        ephem.houses(jd, lat, lon, ord("P"))


def benchmark_time_range(n_days: int = 365) -> None:
    """Benchmark calculations over a time range (simulating ephemeris generation)."""
    start_jd = ephem.julday(2000, 1, 1, 12.0)

    for i in range(n_days):
        jd = start_jd + i
        # Sun and Moon positions
        ephem.calc_ut(jd, SUN, FLG_SPEED)
        ephem.calc_ut(jd, MOON, FLG_SPEED)


def run_all_benchmarks() -> None:
    """Run all benchmarks."""
    print("Running benchmark suite...")
    print("-" * 60)

    benchmark_julian_day(1000)
    benchmark_planet_calculations(100)
    benchmark_house_calculations(100)
    benchmark_lunar_nodes(100)
    benchmark_ayanamsa(100)
    benchmark_full_chart(50)
    benchmark_time_range(365)


def profile_with_cprofile() -> pstats.Stats:
    """Profile benchmarks using cProfile."""
    profiler = cProfile.Profile()
    profiler.enable()

    run_all_benchmarks()

    profiler.disable()
    return pstats.Stats(profiler)


def main():
    """Profile the declared local benchmark workload and print hot paths.

    Iteration counts and cProfile sort orders are engineering diagnostics, not
    scientific parameters.  This command changes no ephemeris data and makes
    no portable performance guarantee for other machines or environments.
    """
    print("=" * 70)
    print("libephemeris Performance Profiling")
    print("=" * 70)
    print()

    # Run profiling
    print("Profiling in progress...")
    start = time.perf_counter()
    stats = profile_with_cprofile()
    elapsed = time.perf_counter() - start
    print(f"Total benchmark time: {elapsed:.2f} seconds")
    print()

    # Show results sorted by cumulative time
    print("=" * 70)
    print("TOP 30 FUNCTIONS BY CUMULATIVE TIME")
    print("=" * 70)
    stream = io.StringIO()
    stats.stream = stream
    stats.sort_stats(SortKey.CUMULATIVE)
    stats.print_stats(30)
    print(stream.getvalue())

    # Show results sorted by total time (self-time)
    print("=" * 70)
    print("TOP 30 FUNCTIONS BY SELF TIME (hot paths)")
    print("=" * 70)
    stream = io.StringIO()
    stats.stream = stream
    stats.sort_stats(SortKey.TIME)
    stats.print_stats(30)
    print(stream.getvalue())

    # Show callers for top functions
    print("=" * 70)
    print("CALL GRAPH FOR KEY FUNCTIONS")
    print("=" * 70)
    stream = io.StringIO()
    stats.stream = stream
    stats.print_callers("_calc_body", "houses", "iau2000b_radians", "at")
    print(stream.getvalue())

    # Summary
    print("=" * 70)
    print("PROFILING SUMMARY")
    print("=" * 70)
    print(f"Total time: {elapsed:.2f} seconds")
    print()
    print("Key insights:")
    print("- Functions with high 'tottime' are candidates for optimization")
    print("- Functions called many times with moderate 'tottime' per call")
    print("  may benefit from caching")
    print("- Functions from Skyfield (at, observe, apparent) are external")
    print("  but their call patterns can be optimized")


if __name__ == "__main__":
    main()
