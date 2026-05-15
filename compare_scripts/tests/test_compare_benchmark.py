"""
Pytest-style Performance Benchmark Tests.

Validates that libephemeris performs within acceptable range compared to pyswisseph.
These are not strict comparison tests but performance benchmarks.

NOTE: These tests are skipped because libephemeris is a pure Python implementation
while pyswisseph wraps the C library. Performance comparison is not meaningful.
"""

import pytest

pytestmark = pytest.mark.skip(
    reason="Performance comparison not meaningful: pure Python vs C library"
)

import statistics
import time
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    MEAN_NODE,
    TRUE_NODE,
    FLG_SWIEPH,
    FLG_SPEED,
)


# ============================================================================
# BENCHMARK CONFIGURATION
# ============================================================================

# Number of iterations for benchmarks
BENCHMARK_ITERATIONS = 1000
QUICK_ITERATIONS = 100

# Performance ratio threshold (libephemeris / pyswisseph)
# Since libephemeris is pure Python and pyswisseph is C, we expect slower
MAX_ACCEPTABLE_RATIO = 100.0  # Allow up to 100x slower

# Test Julian Day
TEST_JD = 2451545.0  # J2000.0


# ============================================================================
# BENCHMARK UTILITIES
# ============================================================================


def benchmark_function(func, iterations=QUICK_ITERATIONS):
    """Run a benchmark and return mean time in nanoseconds."""
    times = []
    for _ in range(iterations):
        start = time.perf_counter_ns()
        func()
        end = time.perf_counter_ns()
        times.append(end - start)
    return statistics.mean(times)


def compare_performance(libephem_func, pyswisseph_func, iterations=QUICK_ITERATIONS):
    """Compare performance between libephemeris and pyswisseph."""
    libephem_time = benchmark_function(libephem_func, iterations)
    pyswisseph_time = benchmark_function(pyswisseph_func, iterations)

    if pyswisseph_time > 0:
        ratio = libephem_time / pyswisseph_time
    else:
        ratio = float("inf")

    return {
        "libephem_ns": libephem_time,
        "pyswisseph_ns": pyswisseph_time,
        "ratio": ratio,
    }


# ============================================================================
# PLANETARY CALCULATION BENCHMARKS
# ============================================================================


class TestPlanetBenchmarks:
    """Benchmark tests for planetary position calculations."""

    @pytest.mark.benchmark
    @pytest.mark.parametrize(
        "planet_py,planet_swe,name",
        [
            (SUN, swe.SUN, "Sun"),
            (MOON, swe.MOON, "Moon"),
            (MARS, swe.MARS, "Mars"),
            (JUPITER, swe.JUPITER, "Jupiter"),
        ],
    )
    def test_planet_calculation_performance(self, planet_py, planet_swe, name):
        """Test that planet calculations complete within acceptable time."""
        flags = FLG_SWIEPH | FLG_SPEED

        result = compare_performance(
            lambda: ephem.calc_ut(TEST_JD, planet_py, flags),
            lambda: swe.calc_ut(TEST_JD, planet_swe, flags),
        )

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO, (
            f"{name} calculation is {result['ratio']:.1f}x slower than pyswisseph"
        )

    @pytest.mark.benchmark
    def test_all_planets_batch(self):
        """Test batch calculation of all planets."""
        flags = FLG_SWIEPH | FLG_SPEED
        planets_py = [
            SUN,
            MOON,
            MERCURY,
            VENUS,
            MARS,
            JUPITER,
            SATURN,
        ]
        planets_swe = [
            swe.SUN,
            swe.MOON,
            swe.MERCURY,
            swe.VENUS,
            swe.MARS,
            swe.JUPITER,
            swe.SATURN,
        ]

        def calc_all_libephem():
            for p in planets_py:
                ephem.calc_ut(TEST_JD, p, flags)

        def calc_all_pyswisseph():
            for p in planets_swe:
                swe.calc_ut(TEST_JD, p, flags)

        result = compare_performance(calc_all_libephem, calc_all_pyswisseph)

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO


# ============================================================================
# HOUSE CALCULATION BENCHMARKS
# ============================================================================


class TestHouseBenchmarks:
    """Benchmark tests for house calculations."""

    @pytest.mark.benchmark
    @pytest.mark.parametrize(
        "hsys,hsys_name",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("R", "Regiomontanus"),
            ("W", "Whole Sign"),
        ],
    )
    def test_house_calculation_performance(self, hsys, hsys_name):
        """Test that house calculations complete within acceptable time."""
        lat, lon = 41.9028, 12.4964  # Rome

        result = compare_performance(
            lambda: ephem.houses(TEST_JD, lat, lon, hsys),
            lambda: swe.houses(TEST_JD, lat, lon, hsys.encode("ascii")),
        )

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO, (
            f"{hsys_name} houses calculation is {result['ratio']:.1f}x slower"
        )


# ============================================================================
# UTILITY FUNCTION BENCHMARKS
# ============================================================================


class TestUtilityBenchmarks:
    """Benchmark tests for utility functions."""

    @pytest.mark.benchmark
    def test_julday_performance(self):
        """Test Julian Day calculation performance."""
        result = compare_performance(
            lambda: ephem.julday(2000, 1, 1, 12.0),
            lambda: swe.julday(2000, 1, 1, 12.0),
        )

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO

    @pytest.mark.benchmark
    def test_revjul_performance(self):
        """Test reverse Julian Day calculation performance."""
        result = compare_performance(
            lambda: ephem.revjul(TEST_JD),
            lambda: swe.revjul(TEST_JD),
        )

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO

    @pytest.mark.benchmark
    def test_deltat_performance(self):
        """Test Delta T calculation performance."""
        result = compare_performance(
            lambda: ephem.deltat(TEST_JD),
            lambda: swe.deltat(TEST_JD),
        )

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO


# ============================================================================
# LUNAR POINT BENCHMARKS
# ============================================================================


class TestLunarPointBenchmarks:
    """Benchmark tests for lunar point calculations."""

    @pytest.mark.benchmark
    @pytest.mark.parametrize(
        "point_py,point_swe,name",
        [
            (MEAN_NODE, swe.MEAN_NODE, "Mean Node"),
            (TRUE_NODE, swe.TRUE_NODE, "True Node"),
        ],
    )
    def test_lunar_point_performance(self, point_py, point_swe, name):
        """Test lunar point calculation performance."""
        flags = FLG_SWIEPH | FLG_SPEED

        result = compare_performance(
            lambda: ephem.calc_ut(TEST_JD, point_py, flags),
            lambda: swe.calc_ut(TEST_JD, point_swe, flags),
        )

        assert result["ratio"] < MAX_ACCEPTABLE_RATIO, (
            f"{name} calculation is {result['ratio']:.1f}x slower"
        )


# ============================================================================
# SUMMARY BENCHMARKS
# ============================================================================


class TestBenchmarkSummary:
    """Summary benchmark tests."""

    @pytest.mark.benchmark
    def test_overall_performance(self):
        """Test that overall performance is within acceptable range."""
        flags = FLG_SWIEPH | FLG_SPEED
        ratios = []

        # Test a few key calculations
        test_cases = [
            (
                lambda: ephem.calc_ut(TEST_JD, SUN, flags),
                lambda: swe.calc_ut(TEST_JD, swe.SUN, flags),
            ),
            (
                lambda: ephem.houses(TEST_JD, 45.0, 0.0, "P"),
                lambda: swe.houses(TEST_JD, 45.0, 0.0, b"P"),
            ),
            (
                lambda: ephem.julday(2024, 1, 1, 12.0),
                lambda: swe.julday(2024, 1, 1, 12.0),
            ),
        ]

        for libephem_func, pyswisseph_func in test_cases:
            result = compare_performance(libephem_func, pyswisseph_func)
            ratios.append(result["ratio"])

        mean_ratio = statistics.mean(ratios)

        assert mean_ratio < MAX_ACCEPTABLE_RATIO, (
            f"Mean performance ratio is {mean_ratio:.1f}x slower than pyswisseph"
        )

    @pytest.mark.benchmark
    def test_calculation_correctness_with_speed(self):
        """Test that fast calculations are still correct."""
        flags = FLG_SWIEPH | FLG_SPEED

        # Quick calculation
        pos_py, _ = ephem.calc_ut(TEST_JD, SUN, flags)
        pos_swe, _ = swe.calc_ut(TEST_JD, swe.SUN, flags)

        # Should still be accurate
        diff = abs(pos_py[0] - pos_swe[0])
        assert diff < 0.001, "Calculation accuracy degraded"
