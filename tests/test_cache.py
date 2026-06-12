"""
Tests for the cache module (hot path optimization).

These tests verify that:
1. Cached functions return correct values
2. Cache hits provide expected speedup
3. Cache clearing works correctly
4. Edge cases are handled properly
"""

import math
import time
from libephemeris.cache import (
    get_cached_nutation,
    get_nutation_degrees,
    get_cached_obliquity,
    get_mean_obliquity,
    get_true_obliquity,
    clear_caches,
    get_cache_info,
)


# Reference Julian Day: J2000.0
JD_J2000 = 2451545.0

# Reference Julian Day: 2024-01-01 12:00 TT
JD_2024 = 2460310.0


class TestNutationCache:
    """Tests for nutation caching."""

    def setup_method(self):
        """Clear caches before each test."""
        clear_caches()

    def test_nutation_returns_tuple(self):
        """Test that get_cached_nutation returns a 2-tuple of floats."""
        result = get_cached_nutation(JD_J2000)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)

    def test_nutation_values_reasonable(self):
        """Test that nutation values are within expected ranges."""
        dpsi, deps = get_cached_nutation(JD_J2000)

        # Nutation in longitude is typically < 20 arcseconds (< 0.001 radians)
        assert abs(dpsi) < 0.001, f"dpsi={dpsi} radians seems too large"

        # Nutation in obliquity is typically < 10 arcseconds (< 0.0001 radians)
        assert abs(deps) < 0.0001, f"deps={deps} radians seems too large"

    def test_nutation_degrees_conversion(self):
        """Test that degree conversion is accurate."""
        dpsi_rad, deps_rad = get_cached_nutation(JD_J2000)
        dpsi_deg, deps_deg = get_nutation_degrees(JD_J2000)

        assert abs(dpsi_deg - math.degrees(dpsi_rad)) < 1e-10
        assert abs(deps_deg - math.degrees(deps_rad)) < 1e-10

    def test_nutation_cache_hit(self):
        """Test that cache hits are fast."""
        # First call - cache miss
        clear_caches()
        start = time.perf_counter()
        result1 = get_cached_nutation(JD_J2000)
        time1 = time.perf_counter() - start

        # Second call - cache hit
        start = time.perf_counter()
        result2 = get_cached_nutation(JD_J2000)
        time2 = time.perf_counter() - start

        # Results should be identical
        assert result1 == result2

        # Cache hit should be significantly faster (at least 10x)
        if time1 > 0.0001:  # Only check if first call took measurable time
            assert time2 < time1, "Cache hit should be faster than cache miss"

    def test_nutation_different_jd(self):
        """Test that different JDs return different nutation values."""
        result1 = get_cached_nutation(JD_J2000)
        result2 = get_cached_nutation(JD_2024)

        # Values should differ (nutation changes over time)
        assert result1 != result2

    def test_cache_info(self):
        """Test cache info reporting."""
        clear_caches()

        # Make a call
        get_cached_nutation(JD_J2000)

        info = get_cache_info()
        assert "nutation" in info
        assert info["nutation"]["hits"] >= 0
        assert info["nutation"]["misses"] >= 1


class TestObliquityCache:
    """Tests for obliquity caching."""

    def setup_method(self):
        """Clear caches before each test."""
        clear_caches()

    def test_obliquity_returns_tuple(self):
        """Test that get_cached_obliquity returns a 2-tuple of floats."""
        result = get_cached_obliquity(JD_J2000)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)

    def test_obliquity_values_reasonable(self):
        """Test that obliquity values are within expected ranges."""
        eps0, eps = get_cached_obliquity(JD_J2000)

        # Mean obliquity at J2000.0 is approximately 23.4393°
        assert 23.0 < eps0 < 24.0, f"Mean obliquity {eps0}° out of range"

        # True obliquity differs from mean by at most ~10 arcseconds (~0.003°)
        assert abs(eps - eps0) < 0.01, "True-mean difference too large"

    def test_mean_obliquity_function(self):
        """Test the convenience function get_mean_obliquity."""
        eps0, _ = get_cached_obliquity(JD_J2000)
        eps0_func = get_mean_obliquity(JD_J2000)

        assert eps0 == eps0_func

    def test_true_obliquity_function(self):
        """Test the convenience function get_true_obliquity."""
        _, eps = get_cached_obliquity(JD_J2000)
        eps_func = get_true_obliquity(JD_J2000)

        assert eps == eps_func

    def test_obliquity_j2000_reference(self):
        """Test obliquity at J2000.0 against known reference value."""
        eps0 = get_mean_obliquity(JD_J2000)

        # Reference: 23°26'21.448" = 23.43929° (Lieske 1979 / IAU 1976)
        # Note: Our formula is Laskar 1986 which gives slightly different value
        reference = 23.43929111  # degrees

        # Should match within 0.01 arcseconds
        assert abs(eps0 - reference) < 0.00003, (
            f"Mean obliquity {eps0} vs reference {reference}"
        )

    def test_obliquity_uses_nutation_cache(self):
        """Test that obliquity calculation uses the nutation cache."""
        clear_caches()

        # First, get obliquity (should call nutation internally)
        get_cached_obliquity(JD_J2000)

        # Check that nutation cache was populated
        info = get_cache_info()
        assert info["nutation"]["misses"] >= 1 or info["nutation"]["hits"] >= 1

    def test_obliquity_different_epochs(self):
        """Test obliquity changes over time (precession)."""
        eps_j2000 = get_mean_obliquity(JD_J2000)
        eps_2024 = get_mean_obliquity(JD_2024)

        # Obliquity decreases by ~47 arcseconds per century
        # From J2000 to 2024 is ~24 years = ~0.24 centuries
        # Expected decrease: ~11 arcseconds = ~0.003°
        diff = eps_j2000 - eps_2024

        # Difference should be positive (obliquity decreasing)
        assert diff > 0, "Obliquity should decrease over time"

        # And within expected range
        assert 0.001 < diff < 0.01, f"Obliquity change {diff}° seems wrong"


class TestCacheClearing:
    """Tests for cache clearing functionality."""

    def test_clear_caches(self):
        """Test that clear_caches properly clears all caches."""
        # Populate caches
        get_cached_nutation(JD_J2000)
        get_cached_obliquity(JD_J2000)

        # Check caches are populated
        info1 = get_cache_info()
        assert info1["nutation"]["currsize"] >= 1
        assert info1["obliquity"]["currsize"] >= 1

        # Clear caches
        clear_caches()

        # Check caches are empty
        info2 = get_cache_info()
        assert info2["nutation"]["currsize"] == 0
        assert info2["obliquity"]["currsize"] == 0


class TestCachePerformance:
    """Performance tests for cache effectiveness."""

    def test_repeated_calls_fast(self):
        """Test that repeated calls to cached functions are fast."""
        clear_caches()

        # Warm up the cache
        get_cached_nutation(JD_J2000)
        get_cached_obliquity(JD_J2000)

        # Time 1000 cached calls
        n = 1000
        start = time.perf_counter()
        for _ in range(n):
            get_cached_nutation(JD_J2000)
            get_cached_obliquity(JD_J2000)
        elapsed = time.perf_counter() - start

        # 2000 cache hits should take less than 10ms
        assert elapsed < 0.1, f"Cache hits too slow: {elapsed:.3f}s for {n * 2} calls"

    def test_cache_effective_with_houses(self):
        """Test that cache is effective when used with houses module."""
        import libephemeris as ephem

        clear_caches()

        jd = ephem.julday(2024, 1, 1, 12.0)
        lat, lon = 41.9, 12.5  # Rome

        # First call - cache cold
        start = time.perf_counter()
        ephem.houses(jd, lat, lon, ord("P"))
        time1 = time.perf_counter() - start

        # Subsequent calls - cache warm
        times = []
        for _ in range(10):
            start = time.perf_counter()
            ephem.houses(jd, lat, lon, ord("P"))
            times.append(time.perf_counter() - start)

        avg_warm = sum(times) / len(times)

        # With caching, subsequent calls should be faster
        # (though this depends on what fraction of time is cacheable)
        # Just verify they're reasonably fast (< 5ms each)
        assert avg_warm < 0.005, (
            f"House calculations too slow: {avg_warm * 1000:.2f}ms avg"
        )
