"""
Tests for compare_benchmark.py performance benchmark module.

Tests verify that:
1. Benchmark functions run without errors
2. Results have expected structure and valid data
3. Statistics calculations are correct
4. Quick benchmarks complete in reasonable time
"""

import sys
from pathlib import Path
import time

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "compare_scripts"))

# Skip entire module if compare_benchmark doesn't exist
pytest.importorskip("compare_benchmark", reason="compare_benchmark module not found")

from compare_benchmark import (
    BenchmarkResult,
    ComparisonResult,
    benchmark_coordinates,
    benchmark_houses,
    benchmark_lunar_points,
    benchmark_planets,
    benchmark_sidereal,
    benchmark_utilities,
    compare_benchmark,
    get_quick_benchmark_stats,
    run_benchmark,
    run_benchmarks,
)


class TestBenchmarkResult:
    """Test the BenchmarkResult dataclass."""

    def test_basic_statistics(self):
        """Test mean, median, min, max calculations."""
        result = BenchmarkResult(
            name="test",
            iterations=5,
            times_ns=[100, 200, 300, 400, 500],
            category="Test",
        )

        assert result.mean_ns == 300.0
        assert result.median_ns == 300.0
        assert result.min_ns == 100
        assert result.max_ns == 500

    def test_stdev_calculation(self):
        """Test standard deviation is calculated correctly."""
        result = BenchmarkResult(
            name="test",
            iterations=5,
            times_ns=[100, 100, 100, 100, 100],
            category="Test",
        )
        # All same values should have 0 stdev
        assert result.stdev_ns == 0.0

    def test_total_time(self):
        """Test total time in milliseconds."""
        result = BenchmarkResult(
            name="test",
            iterations=3,
            times_ns=[1_000_000, 2_000_000, 3_000_000],  # 1ms, 2ms, 3ms
            category="Test",
        )
        assert result.total_ms == 6.0  # 6 milliseconds total

    def test_format_time_nanoseconds(self):
        """Test time formatting for nanoseconds."""
        result = BenchmarkResult(name="test", iterations=1, times_ns=[500])
        assert "ns" in result.format_time(500)

    def test_format_time_microseconds(self):
        """Test time formatting for microseconds."""
        result = BenchmarkResult(name="test", iterations=1, times_ns=[5000])
        assert "us" in result.format_time(5000)

    def test_format_time_milliseconds(self):
        """Test time formatting for milliseconds."""
        result = BenchmarkResult(name="test", iterations=1, times_ns=[5_000_000])
        assert "ms" in result.format_time(5_000_000)

    def test_empty_times(self):
        """Test handling of empty times list."""
        result = BenchmarkResult(
            name="test", iterations=0, times_ns=[], category="Test"
        )
        assert result.mean_ns == 0.0
        assert result.median_ns == 0.0
        assert result.min_ns == 0.0
        assert result.max_ns == 0.0


class TestComparisonResult:
    """Test the ComparisonResult dataclass."""

    def test_ratio_calculation(self):
        """Test ratio is calculated correctly."""
        libephem = BenchmarkResult(
            name="libephem", iterations=1, times_ns=[200], category="Test"
        )
        pyswisseph = BenchmarkResult(
            name="pyswisseph", iterations=1, times_ns=[100], category="Test"
        )

        result = ComparisonResult(
            name="test", libephem=libephem, pyswisseph=pyswisseph, category="Test"
        )

        # libephem is 2x slower
        assert result.ratio == 2.0
        assert "slower" in result.speedup

    def test_ratio_faster(self):
        """Test ratio when libephemeris is faster."""
        libephem = BenchmarkResult(
            name="libephem", iterations=1, times_ns=[50], category="Test"
        )
        pyswisseph = BenchmarkResult(
            name="pyswisseph", iterations=1, times_ns=[100], category="Test"
        )

        result = ComparisonResult(
            name="test", libephem=libephem, pyswisseph=pyswisseph, category="Test"
        )

        assert result.ratio == 0.5
        assert "faster" in result.speedup

    def test_ratio_same_speed(self):
        """Test ratio when speeds are equal."""
        libephem = BenchmarkResult(
            name="libephem", iterations=1, times_ns=[100], category="Test"
        )
        pyswisseph = BenchmarkResult(
            name="pyswisseph", iterations=1, times_ns=[100], category="Test"
        )

        result = ComparisonResult(
            name="test", libephem=libephem, pyswisseph=pyswisseph, category="Test"
        )

        assert result.ratio == 1.0
        assert "same" in result.speedup


class TestRunBenchmark:
    """Test the run_benchmark function."""

    def test_basic_benchmark(self):
        """Test running a basic benchmark."""
        counter = [0]

        def simple_func():
            counter[0] += 1

        result = run_benchmark(
            name="test", func=simple_func, iterations=100, warmup=10, category="Test"
        )

        # Should have run warmup + iterations times
        assert counter[0] == 110
        assert result.iterations == 100
        assert len(result.times_ns) == 100
        assert result.name == "test"
        assert result.category == "Test"

    def test_timing_is_positive(self):
        """Test that all timing measurements are positive."""
        result = run_benchmark(
            name="test",
            func=lambda: sum(range(100)),
            iterations=50,
            warmup=5,
            category="Test",
        )

        assert all(t > 0 for t in result.times_ns)

    def test_benchmark_captures_real_time(self):
        """Test that benchmark captures actual execution time."""
        import time

        def slow_func():
            time.sleep(0.0001)  # 100 microseconds

        result = run_benchmark(
            name="test", func=slow_func, iterations=10, warmup=2, category="Test"
        )

        # Each call should take at least 100 microseconds = 100,000 nanoseconds
        assert result.mean_ns > 50_000  # Allow some tolerance


class TestCompareBenchmark:
    """Test the compare_benchmark function."""

    def test_compare_returns_comparison_result(self):
        """Test that compare_benchmark returns a ComparisonResult."""
        result = compare_benchmark(
            name="test",
            libephem_func=lambda: sum(range(100)),
            pyswisseph_func=lambda: sum(range(100)),
            iterations=50,
            category="Test",
        )

        assert isinstance(result, ComparisonResult)
        assert result.name == "test"
        assert result.category == "Test"
        assert result.libephem is not None
        assert result.pyswisseph is not None
        assert result.ratio > 0


class TestBenchmarkCategories:
    """Test individual benchmark category functions."""

    def test_benchmark_planets(self):
        """Test planet benchmarks run successfully."""
        results = benchmark_planets(iterations=10)

        assert len(results) > 0
        assert all(isinstance(r, ComparisonResult) for r in results)
        assert all(r.category == "Planets" for r in results)

    def test_benchmark_houses(self):
        """Test house benchmarks run successfully."""
        results = benchmark_houses(iterations=10)

        assert len(results) > 0
        assert all(isinstance(r, ComparisonResult) for r in results)
        assert all(r.category == "Houses" for r in results)

    def test_benchmark_lunar_points(self):
        """Test lunar point benchmarks run successfully."""
        results = benchmark_lunar_points(iterations=10)

        assert len(results) > 0
        assert all(isinstance(r, ComparisonResult) for r in results)
        assert all(r.category == "Lunar Points" for r in results)

    def test_benchmark_utilities(self):
        """Test utility benchmarks run successfully."""
        results = benchmark_utilities(iterations=10)

        assert len(results) > 0
        assert all(isinstance(r, ComparisonResult) for r in results)
        assert all(r.category == "Utilities" for r in results)

    def test_benchmark_sidereal(self):
        """Test sidereal benchmarks run successfully."""
        results = benchmark_sidereal(iterations=10)

        assert len(results) > 0
        assert all(isinstance(r, ComparisonResult) for r in results)
        assert all(r.category == "Sidereal" for r in results)

    def test_benchmark_coordinates(self):
        """Test coordinate benchmarks run successfully."""
        results = benchmark_coordinates(iterations=10)

        assert len(results) > 0
        assert all(isinstance(r, ComparisonResult) for r in results)
        assert all(r.category == "Coordinates" for r in results)


class TestRunBenchmarks:
    """Test the run_benchmarks main function."""

    def test_run_all_benchmarks(self):
        """Test running all benchmarks with minimal iterations."""
        result = run_benchmarks(iterations=5, verbose=False)

        assert "results" in result
        assert "summary" in result
        assert "report" in result

        summary = result["summary"]
        assert "total_benchmarks" in summary
        assert "mean_ratio" in summary
        assert "median_ratio" in summary
        assert summary["total_benchmarks"] > 0

    def test_run_single_category(self):
        """Test running a single category of benchmarks."""
        result = run_benchmarks(iterations=5, categories=["utilities"], verbose=False)

        assert len(result["results"]) > 0
        # All results should be from the utilities category
        assert all(r.category == "Utilities" for r in result["results"])

    def test_report_generation(self):
        """Test that a performance report is generated."""
        result = run_benchmarks(iterations=5, categories=["utilities"], verbose=False)

        report = result["report"]
        assert isinstance(report, str)
        assert len(report) > 0
        assert "PERFORMANCE" in report or "Benchmark" in report


class TestGetQuickBenchmarkStats:
    """Test the quick benchmark stats function."""

    def test_quick_stats_returns_summary(self):
        """Test that quick stats returns a summary dict."""
        stats = get_quick_benchmark_stats(iterations=10)

        assert isinstance(stats, dict)
        assert "total_benchmarks" in stats
        assert "mean_ratio" in stats
        assert "median_ratio" in stats
        assert "min_ratio" in stats
        assert "max_ratio" in stats

    def test_quick_stats_valid_values(self):
        """Test that quick stats returns valid values."""
        stats = get_quick_benchmark_stats(iterations=10)

        assert stats["total_benchmarks"] > 0
        assert stats["mean_ratio"] > 0
        assert stats["median_ratio"] > 0
        assert stats["min_ratio"] > 0
        assert stats["max_ratio"] >= stats["min_ratio"]


class TestBenchmarkPerformance:
    """Performance-related tests to ensure benchmarks complete in reasonable time."""

    def test_benchmark_overhead(self):
        """Test that benchmark infrastructure has minimal overhead."""
        start = time.perf_counter()
        result = run_benchmark(
            name="test",
            func=lambda: None,  # Empty function
            iterations=1000,
            warmup=100,
        )
        elapsed = time.perf_counter() - start

        # Should complete in less than 1 second for empty function
        assert elapsed < 1.0

    def test_full_benchmark_suite_completes(self):
        """Test that the full benchmark suite can complete."""
        # Run with minimal iterations to verify everything works
        result = run_benchmarks(iterations=3, verbose=False)

        assert result["summary"]["total_benchmarks"] > 20
        assert result["summary"]["mean_ratio"] > 0
