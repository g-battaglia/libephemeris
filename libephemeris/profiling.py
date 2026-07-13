# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Performance profiling utilities for libephemeris.

This module provides tools to profile hot paths and identify performance
bottlenecks in the library. Uses Python's cProfile for accurate function-level
timing and pstats for analysis.

Main Features:
- Profile planetary position calculations
- Identify the most time-consuming functions
- Generate sorted profiling reports
- Compare performance across different operations

Usage:
    >>> from libephemeris.profiling import profile_planetary_calculations
    >>> report = profile_planetary_calculations(iterations=1000)
    >>> print(report.summary())

    # Or use the context manager for custom profiling:
    >>> from libephemeris.profiling import ProfileContext
    >>> with ProfileContext() as p:
    ...     # your code here
    >>> p.print_stats(top_n=20)
"""

from __future__ import annotations

import cProfile
import io
import pstats
from dataclasses import dataclass, field
from typing import Callable, Optional, Any
from functools import wraps


@dataclass
class FunctionStats:
    """Statistics for a single function from profiling."""

    name: str
    calls: int
    total_time: float  # Total time including subfunctions (cumulative)
    own_time: float  # Time in this function only (not subfunctions)
    callers: int = 0

    @property
    def time_per_call_us(self) -> float:
        """Average time per call in microseconds."""
        if self.calls > 0:
            return (self.own_time / self.calls) * 1_000_000
        return 0.0

    @property
    def cumulative_per_call_us(self) -> float:
        """Average cumulative time per call in microseconds."""
        if self.calls > 0:
            return (self.total_time / self.calls) * 1_000_000
        return 0.0


@dataclass
class ProfileReport:
    """
    Container for profiling results with analysis utilities.

    Attributes:
        stats: pstats.Stats object with raw profiling data
        function_stats: List of FunctionStats sorted by cumulative time
        total_calls: Total number of function calls
        total_time: Total execution time in seconds
        operation_name: Name/description of the profiled operation
    """

    stats: Optional[pstats.Stats] = None
    function_stats: list = field(default_factory=list)
    total_calls: int = 0
    total_time: float = 0.0
    operation_name: str = "Unknown"

    def summary(self, top_n: int = 15) -> str:
        """
        Generate a human-readable summary of the profiling results.

        Args:
            top_n: Number of top functions to include

        Returns:
            Formatted string with profiling summary
        """
        lines = [
            "=" * 80,
            f"PROFILING REPORT: {self.operation_name}",
            "=" * 80,
            f"Total function calls: {self.total_calls:,}",
            f"Total execution time: {self.total_time:.4f} seconds",
            "",
            "TOP FUNCTIONS BY CUMULATIVE TIME:",
            "-" * 80,
            f"{'Function':<50} {'Calls':>8} {'Cumul(us)':>12} {'Own(us)':>12}",
            "-" * 80,
        ]

        for fs in self.function_stats[:top_n]:
            # Truncate long function names
            name = fs.name[:48] + ".." if len(fs.name) > 50 else fs.name
            lines.append(
                f"{name:<50} {fs.calls:>8} "
                f"{fs.cumulative_per_call_us:>12.2f} "
                f"{fs.time_per_call_us:>12.2f}"
            )

        lines.extend(
            [
                "-" * 80,
                "",
                "HOT PATH ANALYSIS:",
                "-" * 40,
            ]
        )

        # Identify key bottlenecks
        hot_paths = self._identify_hot_paths()
        for category, funcs in hot_paths.items():
            if funcs:
                lines.append(f"\n{category}:")
                for func in funcs[:5]:
                    pct = (
                        (func.total_time / self.total_time) * 100
                        if self.total_time > 0
                        else 0
                    )
                    lines.append(f"  - {func.name}: {pct:.1f}% of total time")

        lines.append("=" * 80)
        return "\n".join(lines)

    def _identify_hot_paths(self) -> dict:
        """Group functions by category for hot path analysis."""
        categories: dict[str, list] = {
            "Skyfield Operations": [],
            "Libephemeris Core": [],
            "Math Operations": [],
            "I/O and File": [],
            "Other": [],
        }

        for fs in self.function_stats:
            name_lower = fs.name.lower()
            if any(x in name_lower for x in ["skyfield", "spice", "jpllib", "segment"]):
                categories["Skyfield Operations"].append(fs)
            elif any(
                x in name_lower for x in ["libephemeris", "planets", "lunar", "houses"]
            ):
                categories["Libephemeris Core"].append(fs)
            elif any(x in name_lower for x in ["math", "sin", "cos", "sqrt", "atan"]):
                categories["Math Operations"].append(fs)
            elif any(
                x in name_lower for x in ["read", "write", "open", "load", "file"]
            ):
                categories["I/O and File"].append(fs)
            else:
                categories["Other"].append(fs)

        return categories

    def get_functions_by_name(self, pattern: str) -> list:
        """
        Find functions matching a name pattern.

        Args:
            pattern: Substring to search for in function names

        Returns:
            List of FunctionStats matching the pattern
        """
        return [fs for fs in self.function_stats if pattern.lower() in fs.name.lower()]

    def print_stats(self, top_n: int = 20, sort_by: str = "cumulative") -> None:
        """
        Print detailed profiling statistics.

        Args:
            top_n: Number of top functions to display
            sort_by: Sort key ('cumulative', 'time', 'calls')
        """
        if self.stats:
            self.stats.sort_stats(sort_by)
            self.stats.print_stats(top_n)


class ProfileContext:
    """
    Context manager for profiling code blocks.

    Usage:
        with ProfileContext(name="my operation") as p:
            # code to profile
        report = p.get_report()
        print(report.summary())
    """

    def __init__(self, name: str = "Custom Profile"):
        """
        Initialize profiling context.

        Args:
            name: Name for this profiling session
        """
        self.name = name
        self.profiler = cProfile.Profile()
        self._report: Optional[ProfileReport] = None

    def __enter__(self) -> "ProfileContext":
        """Start profiling."""
        self.profiler.enable()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Stop profiling and generate report."""
        self.profiler.disable()
        self._generate_report()

    def _generate_report(self) -> None:
        """Generate ProfileReport from collected data."""
        # Create stats object
        stream = io.StringIO()
        stats = pstats.Stats(self.profiler, stream=stream)
        stats.sort_stats("cumulative")

        # Extract function statistics
        function_stats = []
        for key, value in stats.stats.items():  # type: ignore[attr-defined]
            filename, line_number, func_name = key
            # value: (ncalls, totcalls, tottime, cumtime, callers)
            ncalls = value[0]
            tottime = value[2]
            cumtime = value[3]

            # Create readable name
            if filename and not filename.startswith("<"):
                # Extract just the filename without full path
                short_file = filename.split("/")[-1]
                name = f"{short_file}:{func_name}"
            else:
                name = func_name

            function_stats.append(
                FunctionStats(
                    name=name,
                    calls=ncalls,
                    total_time=cumtime,
                    own_time=tottime,
                    callers=len(value[4]) if value[4] else 0,
                )
            )

        # Sort by cumulative time
        function_stats.sort(key=lambda x: x.total_time, reverse=True)

        # Calculate totals
        total_calls = sum(fs.calls for fs in function_stats)
        total_time = stats.total_tt if hasattr(stats, "total_tt") else 0.0

        self._report = ProfileReport(
            stats=stats,
            function_stats=function_stats,
            total_calls=total_calls,
            total_time=total_time,
            operation_name=self.name,
        )

    def get_report(self) -> ProfileReport:
        """
        Get the profiling report.

        Returns:
            ProfileReport with analysis results

        Raises:
            RuntimeError: If called before context exits
        """
        if self._report is None:
            raise RuntimeError("Profile not yet complete. Use within 'with' block.")
        return self._report

    def print_stats(self, top_n: int = 20, sort_by: str = "cumulative") -> None:
        """Print stats to stdout."""
        if self._report:
            self._report.print_stats(top_n, sort_by)


def profile_function(func: Callable) -> Callable:
    """
    Decorator to profile a single function.

    Usage:
        @profile_function
        def my_function():
            ...

        # After calling my_function, access:
        report = my_function.last_profile
    """

    @wraps(func)
    def wrapper(*args, **kwargs) -> Any:
        with ProfileContext(name=func.__name__) as p:
            result = func(*args, **kwargs)
        wrapper.last_profile = p.get_report()  # type: ignore[attr-defined]
        return result

    wrapper.last_profile = None  # type: ignore[attr-defined]
    return wrapper


def profile_planetary_calculations(
    iterations: int = 1000,
    include_speed: bool = True,
    planets: Optional[list] = None,
) -> ProfileReport:
    """
    Profile planetary position calculations.

    This runs a realistic workload calculating planetary positions and
    profiles the hot paths.

    Args:
        iterations: Number of date iterations
        include_speed: Include velocity calculations (FLG_SPEED)
        planets: List of planet IDs to calculate (default: all major planets)

    Returns:
        ProfileReport with detailed analysis
    """
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
        FLG_SPEED,
    )

    if planets is None:
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

    flags = FLG_SPEED if include_speed else 0
    jd_start = 2451545.0  # J2000.0

    # Warmup to load ephemeris
    for planet in planets:
        ephem.calc_ut(jd_start, planet, flags)

    # Profile the main calculation loop
    with ProfileContext(
        name=f"Planetary Calculations ({iterations} dates x {len(planets)} planets)"
    ) as p:
        for i in range(iterations):
            jd = jd_start + i
            for planet in planets:
                ephem.calc_ut(jd, planet, flags)

    return p.get_report()


def profile_house_calculations(
    iterations: int = 1000,
    house_systems: Optional[list] = None,
) -> ProfileReport:
    """
    Profile house calculations across different house systems.

    Args:
        iterations: Number of iterations per house system
        house_systems: List of house system codes (default: common systems)

    Returns:
        ProfileReport with analysis
    """
    import libephemeris as ephem

    if house_systems is None:
        house_systems = ["P", "K", "R", "C", "E", "W", "M", "O"]

    jd_start = 2451545.0
    lat, lon = 41.9028, 12.4964  # Rome

    # Warmup
    ephem.houses(jd_start, lat, lon, ord("P"))

    with ProfileContext(
        name=f"House Calculations ({iterations} dates x {len(house_systems)} systems)"
    ) as p:
        for i in range(iterations):
            jd = jd_start + i
            for hsys in house_systems:
                ephem.houses(jd, lat, lon, ord(hsys))

    return p.get_report()


def identify_optimization_opportunities(report: ProfileReport) -> list:
    """
    Analyze a profile report and suggest optimization opportunities.

    Args:
        report: ProfileReport from profiling

    Returns:
        List of optimization suggestions as strings
    """
    suggestions = []

    if not report.function_stats:
        return ["No profiling data available"]

    # Find functions taking more than 10% of total time
    for fs in report.function_stats:
        pct = (fs.total_time / report.total_time) * 100 if report.total_time > 0 else 0
        if pct >= 10:
            if "at" in fs.name.lower() and "skyfield" in fs.name.lower():
                suggestions.append(
                    f"HOTSPOT: {fs.name} ({pct:.1f}% of time) - "
                    "Consider caching ephemeris lookups or reducing observe() calls"
                )
            elif "frame_latlon" in fs.name.lower():
                suggestions.append(
                    f"HOTSPOT: {fs.name} ({pct:.1f}% of time) - "
                    "Coordinate transformations are expensive; consider batching"
                )
            elif "_calc_body" in fs.name.lower():
                suggestions.append(
                    f"HOTSPOT: {fs.name} ({pct:.1f}% of time) - "
                    "Core calculation function; focus optimization here"
                )
            else:
                suggestions.append(
                    f"HOTSPOT: {fs.name} ({pct:.1f}% of time) - "
                    "Potential optimization target"
                )

    # Look for repeated function calls that could be cached
    for fs in report.function_stats[:30]:
        if fs.calls > 1000 and fs.time_per_call_us > 10:
            suggestions.append(
                f"CACHE CANDIDATE: {fs.name} called {fs.calls} times, "
                f"{fs.time_per_call_us:.1f}us per call"
            )

    # Check for I/O in hot path
    io_funcs = report.get_functions_by_name("read")
    io_funcs.extend(report.get_functions_by_name("load"))
    for fs in io_funcs:
        if fs.calls > 10:
            suggestions.append(
                f"I/O IN LOOP: {fs.name} called {fs.calls} times - "
                "Ensure data is cached"
            )

    if not suggestions:
        suggestions.append("No obvious optimization opportunities identified")

    return suggestions


def run_full_profile(iterations: int = 500, verbose: bool = True) -> dict:
    """
    Run comprehensive profiling across all major operations.

    Args:
        iterations: Number of iterations for each operation
        verbose: Print results to stdout

    Returns:
        Dictionary with reports for each operation category
    """
    reports = {}

    # Profile planetary calculations
    if verbose:
        print("Profiling planetary calculations...")
    reports["planets"] = profile_planetary_calculations(iterations)

    # Profile house calculations
    if verbose:
        print("Profiling house calculations...")
    reports["houses"] = profile_house_calculations(iterations)

    if verbose:
        print("\n" + "=" * 80)
        print("COMPREHENSIVE PROFILING RESULTS")
        print("=" * 80)

        for name, report in reports.items():
            print(f"\n{report.summary()}")

            suggestions = identify_optimization_opportunities(report)
            if (
                suggestions
            ):  # pragma: no branch - identify_* always returns a non-empty list
                print("\nOPTIMIZATION SUGGESTIONS:")
                for s in suggestions:
                    print(f"  - {s}")

    return reports
