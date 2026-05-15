"""
Tests for the profiling module.

These tests verify that the profiling utilities work correctly
and produce meaningful reports for performance analysis.
"""

import pytest
from libephemeris.profiling import (
    ProfileContext,
    ProfileReport,
    FunctionStats,
    profile_planetary_calculations,
    profile_house_calculations,
    identify_optimization_opportunities,
    profile_function,
)


class TestFunctionStats:
    """Tests for FunctionStats dataclass."""

    def test_time_per_call_calculation(self):
        """Test time per call calculation."""
        stats = FunctionStats(
            name="test_func",
            calls=100,
            total_time=1.0,  # 1 second cumulative
            own_time=0.5,  # 0.5 seconds own time
        )
        # own_time / calls * 1_000_000 = 0.5 / 100 * 1_000_000 = 5000 us
        assert stats.time_per_call_us == 5000.0

    def test_cumulative_per_call_calculation(self):
        """Test cumulative time per call calculation."""
        stats = FunctionStats(
            name="test_func",
            calls=100,
            total_time=1.0,  # 1 second
            own_time=0.5,
        )
        # total_time / calls * 1_000_000 = 1.0 / 100 * 1_000_000 = 10000 us
        assert stats.cumulative_per_call_us == 10000.0

    def test_zero_calls_handling(self):
        """Test handling of zero calls."""
        stats = FunctionStats(
            name="test_func",
            calls=0,
            total_time=0.0,
            own_time=0.0,
        )
        assert stats.time_per_call_us == 0.0
        assert stats.cumulative_per_call_us == 0.0


class TestProfileContext:
    """Tests for ProfileContext context manager."""

    def test_basic_profiling(self):
        """Test basic profiling of a simple function."""
        with ProfileContext(name="test") as p:
            # Simple computation to profile
            total = sum(range(1000))

        report = p.get_report()
        assert report is not None
        assert report.operation_name == "test"
        assert report.total_calls > 0
        assert report.total_time >= 0

    def test_profiling_with_function_calls(self):
        """Test profiling captures function calls."""

        def inner_function():
            return sum(range(100))

        with ProfileContext(name="with_calls") as p:
            for _ in range(10):
                inner_function()

        report = p.get_report()
        # Should capture the inner function calls
        assert report.total_calls >= 10

    def test_report_before_exit_raises(self):
        """Test that getting report before context exit raises."""
        p = ProfileContext(name="test")
        with pytest.raises(RuntimeError):
            p.get_report()

    def test_summary_generation(self):
        """Test that summary is generated correctly."""
        with ProfileContext(name="summary_test") as p:
            _ = [x * 2 for x in range(1000)]

        report = p.get_report()
        summary = report.summary()

        assert "PROFILING REPORT" in summary
        assert "summary_test" in summary
        assert "Total function calls:" in summary
        assert "Total execution time:" in summary

    def test_function_stats_populated(self):
        """Test that function statistics are populated."""
        with ProfileContext(name="stats_test") as p:
            for i in range(100):
                _ = i**2

        report = p.get_report()
        assert len(report.function_stats) > 0

        # Check that FunctionStats have required fields
        for fs in report.function_stats[:5]:
            assert fs.name is not None
            assert isinstance(fs.calls, int)
            assert fs.calls >= 0
            assert isinstance(fs.total_time, float)


class TestProfileReport:
    """Tests for ProfileReport functionality."""

    def test_get_functions_by_name(self):
        """Test filtering functions by name pattern."""
        with ProfileContext(name="filter_test") as p:
            result = sum(range(1000))

        report = p.get_report()
        # Look for built-in functions
        builtins = report.get_functions_by_name("built")
        assert isinstance(builtins, list)

    def test_hot_path_categorization(self):
        """Test that hot paths are categorized."""
        with ProfileContext(name="category_test") as p:
            import math

            for i in range(100):
                _ = math.sin(i)

        report = p.get_report()
        categories = report._identify_hot_paths()

        # Should have category dictionaries
        assert "Skyfield Operations" in categories
        assert "Libephemeris Core" in categories
        assert "Math Operations" in categories
        assert "Other" in categories


class TestProfileFunction:
    """Tests for the profile_function decorator."""

    def test_decorator_preserves_function(self):
        """Test that decorator preserves function behavior."""

        @profile_function
        def add_numbers(a, b):
            return a + b

        result = add_numbers(2, 3)
        assert result == 5

    def test_decorator_creates_profile(self):
        """Test that decorator creates a profile report."""

        @profile_function
        def compute_sum():
            return sum(range(100))

        compute_sum()

        assert compute_sum.last_profile is not None
        assert isinstance(compute_sum.last_profile, ProfileReport)


class TestPlanetaryProfiling:
    """Integration tests for planetary calculation profiling."""

    @pytest.mark.slow
    def test_profile_planetary_calculations(self):
        """Test profiling of planetary calculations."""
        report = profile_planetary_calculations(iterations=10)

        assert report is not None
        assert report.total_time > 0
        assert report.total_calls > 0

        # Should have captured calc_ut calls
        calc_funcs = report.get_functions_by_name("calc")
        assert len(calc_funcs) > 0

    @pytest.mark.slow
    def test_profile_planetary_with_planets_subset(self):
        """Test profiling with specific planets."""
        from libephemeris.constants import SUN, MOON

        report = profile_planetary_calculations(
            iterations=10,
            planets=[SUN, MOON],
        )

        assert report is not None
        assert report.total_time > 0


class TestHouseProfiling:
    """Integration tests for house calculation profiling."""

    @pytest.mark.slow
    def test_profile_house_calculations(self):
        """Test profiling of house calculations."""
        report = profile_house_calculations(iterations=10)

        assert report is not None
        assert report.total_time > 0
        assert report.total_calls > 0

    @pytest.mark.slow
    def test_profile_house_with_systems_subset(self):
        """Test profiling with specific house systems."""
        report = profile_house_calculations(
            iterations=10,
            house_systems=["P", "K"],
        )

        assert report is not None
        assert report.total_time > 0


class TestOptimizationSuggestions:
    """Tests for optimization suggestion generation."""

    def test_suggestions_returned(self):
        """Test that suggestions are generated from report."""
        # Create a mock-ish report with known stats
        report = ProfileReport(
            total_time=10.0,
            function_stats=[
                FunctionStats(
                    name="hot_function",
                    calls=1000,
                    total_time=5.0,  # 50% of total
                    own_time=5.0,
                ),
                FunctionStats(
                    name="skyfield.jpllib._at",
                    calls=500,
                    total_time=3.0,  # 30% of total
                    own_time=3.0,
                ),
            ],
        )

        suggestions = identify_optimization_opportunities(report)
        assert len(suggestions) > 0

    def test_empty_report_suggestions(self):
        """Test suggestions for empty report."""
        report = ProfileReport(total_time=0.0, function_stats=[])

        suggestions = identify_optimization_opportunities(report)
        assert len(suggestions) == 1
        assert "No profiling data" in suggestions[0]


class TestProfileReportSummary:
    """Tests for profile report summary generation."""

    def test_summary_format(self):
        """Test summary output format."""
        report = ProfileReport(
            operation_name="Test Operation",
            total_calls=1000,
            total_time=1.5,
            function_stats=[
                FunctionStats(
                    name="func1",
                    calls=100,
                    total_time=0.5,
                    own_time=0.3,
                ),
                FunctionStats(
                    name="func2",
                    calls=50,
                    total_time=0.3,
                    own_time=0.2,
                ),
            ],
        )

        summary = report.summary(top_n=5)

        # Check header elements
        assert "Test Operation" in summary
        assert "1,000" in summary or "1000" in summary
        assert "1.5" in summary

        # Check function table
        assert "func1" in summary
        assert "func2" in summary
