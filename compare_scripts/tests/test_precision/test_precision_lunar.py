"""
Tests for the compare_lunar.py comparison script.

Validates that the True Node precision comparison function works correctly.
"""

import sys
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestCompareLunarScript:
    """Tests for the compare_lunar comparison script functions."""

    def test_compare_true_node_precision_runs(self):
        """Test that compare_true_node_precision runs without errors."""
        from compare_lunar import compare_true_node_precision

        # Run the comparison function
        stats = compare_true_node_precision()

        # Verify we got results
        assert stats.total > 0, "Should have tested some dates"
        assert stats.passed > 0, "Should have some passing tests"
        assert stats.errors == 0, "Should not have any errors"

    def test_compare_true_node_precision_pass_rate(self):
        """Test that compare_true_node_precision achieves acceptable pass rate."""
        from compare_lunar import compare_true_node_precision

        stats = compare_true_node_precision()

        # Should pass at least 80% of tests
        pass_rate = stats.pass_rate()
        assert pass_rate >= 80.0, f"Pass rate {pass_rate:.1f}% should be >= 80%"

    def test_compare_true_node_precision_test_count(self):
        """Test that compare_true_node_precision tests ~1000 dates."""
        from compare_lunar import compare_true_node_precision

        stats = compare_true_node_precision()

        # Should test approximately 1000 dates (allowing for skipped dates)
        # Some dates are skipped due to ephemeris range, so we allow a range
        assert stats.total >= 700, f"Should test at least 700 dates, got {stats.total}"
        assert stats.total <= 1000, f"Should test at most 1000 dates, got {stats.total}"

    def test_generate_random_jd(self):
        """Test the random Julian Day generator function."""
        from compare_lunar import generate_random_jd

        dates = generate_random_jd(1900, 2100, count=100, seed=42)

        # Verify correct number of dates
        assert len(dates) == 100, f"Should generate 100 dates, got {len(dates)}"

        # Verify tuple structure
        for year, month, day, hour, jd in dates:
            assert 1900 <= year <= 2100, f"Year {year} out of range"
            assert 1 <= month <= 12, f"Month {month} out of range"
            assert 1 <= day <= 28, f"Day {day} out of range"
            assert 0 <= hour < 24, f"Hour {hour} out of range"
            assert jd > 0, f"Julian Day {jd} should be positive"

    def test_generate_random_jd_reproducible(self):
        """Test that random Julian Day generation is reproducible with seed."""
        from compare_lunar import generate_random_jd

        dates1 = generate_random_jd(1900, 2100, count=10, seed=123)
        dates2 = generate_random_jd(1900, 2100, count=10, seed=123)

        assert dates1 == dates2, "Same seed should produce same dates"

    def test_compare_lunar_nodes(self):
        """Test the compare_lunar_nodes function."""
        import swisseph as swe
        from compare_lunar import compare_lunar_nodes

        jd = swe.julday(2000, 1, 1, 12.0)
        results = compare_lunar_nodes(jd, "Test", "2000-01-01")

        assert "mean_node" in results, "Should have mean_node result"
        assert "true_node" in results, "Should have true_node result"

        # Both should have (passed, diff) tuple
        for key, (passed, diff) in results.items():
            # passed can be bool or numpy bool
            assert bool(passed) in (True, False), f"{key} passed should be boolean-like"
            assert isinstance(diff, (float, int)), f"{key} diff should be numeric"
            assert diff >= 0, f"{key} diff should be non-negative"

    def test_compare_lilith(self):
        """Test the compare_lilith function."""
        import swisseph as swe
        from compare_lunar import compare_lilith

        jd = swe.julday(2000, 1, 1, 12.0)
        results = compare_lilith(jd, "Test", "2000-01-01")

        assert "mean_lilith" in results, "Should have mean_lilith result"
        assert "true_lilith" in results, "Should have true_lilith result"

        # Both should have (passed, diff) tuple
        for key, (passed, diff) in results.items():
            # passed can be bool or numpy bool
            assert bool(passed) in (True, False), f"{key} passed should be boolean-like"
            assert isinstance(diff, (float, int)), f"{key} diff should be numeric"
            assert diff >= 0, f"{key} diff should be non-negative"

    def test_main_returns_success(self):
        """Test that main() returns 0 for success."""
        from compare_lunar import main

        result = main()

        # Main should return 0 if pass rate >= 80%
        assert result == 0, f"main() should return 0 (success), got {result}"
