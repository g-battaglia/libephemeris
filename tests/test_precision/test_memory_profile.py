"""
Memory profiling tests for libephemeris.

This module profiles memory usage during long-running calculations to verify:
1. No memory leaks during extended calculation loops (e.g., 10000 positions)
2. Skyfield ephemeris cache is managed correctly
3. Singleton patterns (Loader, SpiceKernel, Timescale) don't leak memory
4. Angles cache is properly managed

Uses Python's tracemalloc for memory tracking.
"""

import gc
import tracemalloc
from dataclasses import dataclass

import pytest

import libephemeris as ephem
from libephemeris import close as ephem_close
from libephemeris import EphemerisContext
from libephemeris.state import (
    get_planets,
    get_timescale,
    get_loader,
    set_angles_cache,
    get_angles_cache,
    clear_angles_cache,
)
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


# Memory growth threshold - maximum allowed memory increase during loops
# A reasonable threshold for 10000 calculations is ~50MB
MAX_MEMORY_GROWTH_MB = 50.0

# Iterations for memory profiling tests
LOOP_ITERATIONS = 10000


@dataclass
class MemoryStats:
    """Memory statistics from a profiling run."""

    peak_mb: float
    current_mb: float
    start_mb: float
    growth_mb: float
    traced_objects: int

    def __str__(self) -> str:
        return (
            f"Memory: start={self.start_mb:.2f}MB, "
            f"current={self.current_mb:.2f}MB, "
            f"peak={self.peak_mb:.2f}MB, "
            f"growth={self.growth_mb:.2f}MB"
        )


def get_memory_stats(start_size: int = 0) -> MemoryStats:
    """Get current memory statistics from tracemalloc."""
    current, peak = tracemalloc.get_traced_memory()
    start_mb = start_size / (1024 * 1024)
    current_mb = current / (1024 * 1024)
    peak_mb = peak / (1024 * 1024)

    # Get number of traced objects
    snapshot = tracemalloc.take_snapshot()
    traced_objects = len(snapshot.statistics("lineno"))

    return MemoryStats(
        peak_mb=peak_mb,
        current_mb=current_mb,
        start_mb=start_mb,
        growth_mb=current_mb - start_mb,
        traced_objects=traced_objects,
    )


@pytest.fixture
def memory_profiler():
    """
    Fixture to set up and tear down memory profiling.

    Usage:
        def test_memory(memory_profiler):
            start_size = memory_profiler.start()
            # ... do work ...
            stats = memory_profiler.stop()
            assert stats.growth_mb < 50
    """

    class MemoryProfiler:
        def __init__(self):
            self._start_size: int = 0
            self._running: bool = False

        def start(self) -> int:
            """Start memory profiling and return initial memory usage."""
            # Force garbage collection before starting
            gc.collect()
            gc.collect()
            gc.collect()

            # Stop any previous tracing
            if tracemalloc.is_tracing():
                tracemalloc.stop()

            # Start fresh tracing
            tracemalloc.start()
            self._start_size, _ = tracemalloc.get_traced_memory()
            self._running = True
            return self._start_size

        def stop(self) -> MemoryStats:
            """Stop profiling and return memory statistics."""
            if not self._running:
                raise RuntimeError("Profiler not started")

            # Force GC before measuring final memory
            gc.collect()
            gc.collect()
            gc.collect()

            stats = get_memory_stats(self._start_size)
            tracemalloc.stop()
            self._running = False
            return stats

        def checkpoint(self) -> MemoryStats:
            """Take a memory snapshot without stopping profiling."""
            if not self._running:
                raise RuntimeError("Profiler not started")
            gc.collect()
            return get_memory_stats(self._start_size)

    return MemoryProfiler()


@pytest.fixture
def ensure_clean_state():
    """Ensure ephemeris state is clean before and after tests."""
    # Clean up before test
    ephem_close()
    gc.collect()

    yield

    # Clean up after test
    ephem_close()
    gc.collect()


class TestMemoryLongCalculations:
    """
    Memory profiling for long calculation loops.

    These tests verify that running many calculations in a loop
    does not cause memory to grow unboundedly.
    """

    @pytest.mark.slow
    def test_calc_ut_10000_positions_no_memory_leak(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Test that 10000 planetary position calculations don't leak memory.

        This test calculates positions for all major planets over 1000 dates,
        totaling 10000 calculations, and verifies memory doesn't grow excessively.
        """
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
        jd_start = 2451545.0  # J2000.0
        dates = [jd_start + i for i in range(1000)]

        progress = progress_reporter("Memory profiling calc_ut", len(dates))

        # Warm up - load ephemeris and timescale
        _ = ephem.calc_ut(jd_start, SUN, 0)
        gc.collect()

        # Start memory profiling
        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(f"MEMORY PROFILE: calc_ut - {len(planets) * len(dates)} calculations")
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        # Run calculations
        for i, jd in enumerate(dates):
            for planet in planets:
                _ = ephem.calc_ut(jd, planet, 0)

            # Report progress at intervals
            if (i + 1) % 100 == 0:
                progress.update(i, f"jd={jd:.1f}")

        progress.done()

        # Get final memory stats
        final_stats = memory_profiler.stop()

        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")

        # Assert memory growth is within acceptable limits
        assert final_stats.growth_mb < MAX_MEMORY_GROWTH_MB, (
            f"Memory grew by {final_stats.growth_mb:.2f}MB during 10000 calculations, "
            f"exceeds limit of {MAX_MEMORY_GROWTH_MB}MB"
        )

    @pytest.mark.slow
    def test_calc_ut_with_speed_memory(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Test memory usage for calc_ut with FLG_SPEED flag.

        Speed calculations involve additional velocity computations
        which may have different memory characteristics.
        """
        planets = [SUN, MOON, MERCURY, MARS, JUPITER]
        jd_start = 2451545.0
        dates = [jd_start + i for i in range(2000)]

        progress = progress_reporter("Memory profiling calc_ut+SPEED", len(dates))

        # Warm up
        _ = ephem.calc_ut(jd_start, SUN, FLG_SPEED)
        gc.collect()

        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(
            f"MEMORY PROFILE: calc_ut + FLG_SPEED - "
            f"{len(planets) * len(dates)} calculations"
        )
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        for i, jd in enumerate(dates):
            for planet in planets:
                _ = ephem.calc_ut(jd, planet, FLG_SPEED)

            if (i + 1) % 200 == 0:
                progress.update(i)

        progress.done()

        final_stats = memory_profiler.stop()
        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")

        assert final_stats.growth_mb < MAX_MEMORY_GROWTH_MB, (
            f"Memory grew by {final_stats.growth_mb:.2f}MB, "
            f"exceeds limit of {MAX_MEMORY_GROWTH_MB}MB"
        )

    @pytest.mark.slow
    def test_houses_10000_calculations_memory(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Test memory usage for 10000 house calculations.

        House calculations involve trigonometric operations and
        may have different memory characteristics than planet calculations.
        """
        house_systems = [ord("P"), ord("K"), ord("E"), ord("W"), ord("R")]
        locations = [
            (41.9, 12.5),  # Rome
            (51.5, -0.1),  # London
            (40.7, -74.0),  # New York
            (-33.9, 151.2),  # Sydney
        ]
        jd_start = 2451545.0
        dates = [jd_start + i for i in range(500)]

        total_calcs = len(house_systems) * len(locations) * len(dates)
        progress = progress_reporter("Memory profiling houses", len(dates))

        # Warm up
        _ = ephem.houses(jd_start, 41.9, 12.5, ord("P"))
        gc.collect()

        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(f"MEMORY PROFILE: houses - {total_calcs} calculations")
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        for i, jd in enumerate(dates):
            for lat, lon in locations:
                for hsys in house_systems:
                    _ = ephem.houses(jd, lat, lon, hsys)

            if (i + 1) % 50 == 0:
                progress.update(i)

        progress.done()

        final_stats = memory_profiler.stop()
        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")

        assert final_stats.growth_mb < MAX_MEMORY_GROWTH_MB, (
            f"Memory grew by {final_stats.growth_mb:.2f}MB, "
            f"exceeds limit of {MAX_MEMORY_GROWTH_MB}MB"
        )


class TestEphemerisCacheManagement:
    """
    Tests for Skyfield ephemeris cache management.

    Verifies that ephemeris files are loaded once and shared correctly,
    and that the singleton pattern doesn't cause memory issues.
    """

    def test_ephemeris_singleton_loads_once(self, ensure_clean_state):
        """
        Test that ephemeris is loaded only once (singleton pattern).

        Multiple calls to get_planets() should return the same object.
        """
        # First call - loads ephemeris
        planets1 = get_planets()
        planets1_id = id(planets1)

        # Second call - should return cached instance
        planets2 = get_planets()
        planets2_id = id(planets2)

        assert planets1_id == planets2_id, (
            "get_planets() should return the same instance (singleton pattern)"
        )

    def test_timescale_singleton_loads_once(self, ensure_clean_state):
        """
        Test that timescale is loaded only once (singleton pattern).
        """
        ts1 = get_timescale()
        ts1_id = id(ts1)

        ts2 = get_timescale()
        ts2_id = id(ts2)

        assert ts1_id == ts2_id, (
            "get_timescale() should return the same instance (singleton pattern)"
        )

    def test_loader_singleton_loads_once(self, ensure_clean_state):
        """
        Test that loader is loaded only once (singleton pattern).
        """
        loader1 = get_loader()
        loader1_id = id(loader1)

        loader2 = get_loader()
        loader2_id = id(loader2)

        assert loader1_id == loader2_id, (
            "get_loader() should return the same instance (singleton pattern)"
        )

    def test_close_releases_resources(self, ensure_clean_state, memory_profiler):
        """
        Test that close() properly releases ephemeris resources.

        After close(), the next call should create new instances.
        """
        # Load ephemeris
        planets1 = get_planets()
        planets1_id = id(planets1)

        # Close and release
        ephem_close()
        gc.collect()

        # Load again - should be new instance
        planets2 = get_planets()
        planets2_id = id(planets2)

        assert planets1_id != planets2_id, (
            "After close(), get_planets() should return a new instance"
        )

    @pytest.mark.slow
    def test_repeated_close_and_reload_no_leak(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Test that repeatedly closing and reloading ephemeris doesn't leak memory.

        This tests a common pattern in long-running applications.
        """
        iterations = 50
        progress = progress_reporter("Close/reload cycles", iterations)

        # Initial load and close
        _ = get_planets()
        ephem_close()
        gc.collect()

        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(f"MEMORY PROFILE: Close/Reload cycles - {iterations} iterations")
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        for i in range(iterations):
            # Load ephemeris
            planets = get_planets()
            ts = get_timescale()

            # Do some calculations
            for jd in range(2451545, 2451555):
                _ = ephem.calc_ut(jd, SUN, 0)

            # Close and cleanup
            ephem_close()
            gc.collect()

            progress.update(i)

        progress.done()

        # Force final GC
        gc.collect()
        gc.collect()

        final_stats = memory_profiler.stop()
        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")

        # Allow more memory for repeated load/unload cycles
        max_growth = MAX_MEMORY_GROWTH_MB * 2
        assert final_stats.growth_mb < max_growth, (
            f"Memory grew by {final_stats.growth_mb:.2f}MB during {iterations} "
            f"close/reload cycles, exceeds limit of {max_growth}MB"
        )


class TestAnglesCacheManagement:
    """
    Tests for angles cache management.

    Verifies that the angles cache used for Arabic parts calculations
    is properly managed and doesn't leak memory.
    """

    def test_angles_cache_set_and_get(self):
        """Test basic set/get operations on angles cache."""
        # Clear any existing cache
        clear_angles_cache()

        # Set cache
        test_angles = {"Sun": 120.5, "Moon": 240.3, "Asc": 15.7}
        set_angles_cache(test_angles)

        # Get cache
        cached = get_angles_cache()

        assert cached == test_angles, "Cached angles should match set values"

    def test_angles_cache_clear(self):
        """Test that clear_angles_cache() properly empties the cache."""
        # Set some values
        set_angles_cache({"Sun": 100.0, "Moon": 200.0})

        # Clear
        clear_angles_cache()

        # Verify empty
        cached = get_angles_cache()
        assert len(cached) == 0, "Cache should be empty after clear"

    def test_angles_cache_set_makes_copy(self):
        """Test that set_angles_cache creates a copy to prevent mutation."""
        original = {"Sun": 120.5, "Moon": 240.3}
        set_angles_cache(original)

        # Mutate original
        original["Sun"] = 999.9

        # Cache should be unaffected
        cached = get_angles_cache()
        assert cached["Sun"] == 120.5, "Cache should not be affected by mutation"

    @pytest.mark.slow
    def test_repeated_cache_updates_no_leak(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Test that repeatedly updating angles cache doesn't leak memory.

        This simulates calculating many charts where the cache is updated
        for each chart.
        """
        iterations = 10000
        progress = progress_reporter("Angles cache updates", iterations)

        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(f"MEMORY PROFILE: Angles cache updates - {iterations} iterations")
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        for i in range(iterations):
            # Create new angles dict (simulating chart calculation)
            angles = {
                "Sun": float(i % 360),
                "Moon": float((i * 13) % 360),
                "Mercury": float((i * 5) % 360),
                "Venus": float((i * 7) % 360),
                "Mars": float((i * 11) % 360),
                "Asc": float((i * 2) % 360),
                "MC": float((i * 3) % 360),
            }
            set_angles_cache(angles)

            # Clear after use (typical pattern)
            if i % 10 == 0:
                clear_angles_cache()

            if (i + 1) % 1000 == 0:
                progress.update(i)

        progress.done()

        # Final clear and GC
        clear_angles_cache()
        gc.collect()

        final_stats = memory_profiler.stop()
        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")

        # Cache operations should have minimal memory footprint
        assert final_stats.growth_mb < 10.0, (
            f"Memory grew by {final_stats.growth_mb:.2f}MB during cache updates, "
            "cache may be leaking"
        )


class TestEphemerisContextMemory:
    """
    Tests for EphemerisContext memory management.

    Verifies that the thread-safe EphemerisContext properly shares
    resources and doesn't leak memory.
    """

    def test_context_shares_ephemeris(self, ensure_clean_state):
        """
        Test that multiple EphemerisContext instances share ephemeris data.
        """
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        planets1 = ctx1.get_planets()
        planets2 = ctx2.get_planets()

        assert id(planets1) == id(planets2), (
            "Multiple contexts should share the same ephemeris data"
        )

    def test_context_shares_timescale(self, ensure_clean_state):
        """
        Test that multiple EphemerisContext instances share timescale.
        """
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ts1 = ctx1.get_timescale()
        ts2 = ctx2.get_timescale()

        assert id(ts1) == id(ts2), "Multiple contexts should share the same timescale"

    def test_context_has_independent_state(self, ensure_clean_state):
        """
        Test that EphemerisContext instances have independent state.
        """
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        # Set different locations
        ctx1.set_topo(12.5, 41.9, 0)  # Rome
        ctx2.set_topo(-74.0, 40.7, 0)  # New York

        # Verify independence
        topo1 = ctx1.get_topo()
        topo2 = ctx2.get_topo()

        assert topo1 is not topo2, "Contexts should have independent topo state"

    @pytest.mark.slow
    def test_many_contexts_no_leak(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Test that creating many EphemerisContext instances doesn't leak memory.

        This simulates a multi-threaded application creating contexts per request.
        """
        iterations = 1000
        progress = progress_reporter("Creating contexts", iterations)

        # Initial context to load shared resources
        _ = EphemerisContext()
        gc.collect()

        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(f"MEMORY PROFILE: EphemerisContext creation - {iterations} instances")
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        contexts = []
        for i in range(iterations):
            ctx = EphemerisContext()
            # Use valid coordinate ranges: lon [-180, 180], lat [-90, 90]
            lon = float(i % 360 - 180)  # -180 to +179
            lat = float(i % 180 - 90)  # -90 to +89
            ctx.set_topo(lon, lat, 0)
            ctx.set_sid_mode(i % 43)  # Cycle through ayanamshas

            # Do some calculations
            _ = ctx.calc_ut(2451545.0 + i, SUN, 0)

            # Keep reference for some, let others be GC'd
            if i % 10 == 0:
                contexts.append(ctx)

            if (i + 1) % 100 == 0:
                gc.collect()
                progress.update(i)

        progress.done()

        # Clear references and collect
        contexts.clear()
        gc.collect()
        gc.collect()

        final_stats = memory_profiler.stop()
        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")

        # Contexts should be GC'd properly
        assert final_stats.growth_mb < MAX_MEMORY_GROWTH_MB, (
            f"Memory grew by {final_stats.growth_mb:.2f}MB during context creation, "
            f"exceeds limit of {MAX_MEMORY_GROWTH_MB}MB"
        )


class TestMemorySummary:
    """
    Summary test that runs a comprehensive memory profile.
    """

    @pytest.mark.slow
    def test_comprehensive_memory_profile(
        self, memory_profiler, ensure_clean_state, progress_reporter
    ):
        """
        Comprehensive memory profile simulating real-world usage.

        This test simulates calculating many charts with all typical operations:
        - Planet positions
        - House cusps
        - Ayanamsha
        - Cache operations
        """
        num_charts = 500
        progress = progress_reporter("Comprehensive memory test", num_charts)

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

        # Warm up
        jd = 2451545.0
        for planet in planets:
            _ = ephem.calc_ut(jd, planet, FLG_SPEED)
        _ = ephem.houses(jd, 41.9, 12.5, ord("P"))
        gc.collect()

        memory_profiler.start()
        initial_stats = memory_profiler.checkpoint()

        print(f"\n{'=' * 70}")
        print(
            f"COMPREHENSIVE MEMORY PROFILE: Simulating {num_charts} chart calculations"
        )
        print(f"{'=' * 70}")
        print(f"Initial: {initial_stats}")

        for i in range(num_charts):
            jd = 2451545.0 + i * 10  # Different dates

            # Calculate all planets with speed
            for planet in planets:
                _ = ephem.calc_ut(jd, planet, FLG_SPEED)

            # Calculate houses for different locations
            for lat, lon in [(41.9, 12.5), (51.5, -0.1), (40.7, -74.0)]:
                _ = ephem.houses(jd, lat, lon, ord("P"))

            # Update angles cache
            angles = {
                "Sun": float(i % 360),
                "Moon": float((i * 13) % 360),
                "Asc": float((i * 2) % 360),
            }
            set_angles_cache(angles)

            # Clear cache periodically
            if i % 50 == 0:
                clear_angles_cache()
                gc.collect()

            if (i + 1) % 50 == 0:
                checkpoint = memory_profiler.checkpoint()
                progress.update(i, f"mem={checkpoint.current_mb:.1f}MB")

        progress.done()

        clear_angles_cache()
        gc.collect()
        gc.collect()

        final_stats = memory_profiler.stop()

        print(f"Final: {final_stats}")
        print(f"{'=' * 70}")
        print("\nSUMMARY:")
        print(f"  Charts calculated: {num_charts}")
        print(f"  Total calculations: {num_charts * (len(planets) + 3)} per chart")
        print(f"  Memory start: {final_stats.start_mb:.2f}MB")
        print(f"  Memory end: {final_stats.current_mb:.2f}MB")
        print(f"  Memory growth: {final_stats.growth_mb:.2f}MB")
        print(f"  Peak memory: {final_stats.peak_mb:.2f}MB")
        print(f"{'=' * 70}")

        assert final_stats.growth_mb < MAX_MEMORY_GROWTH_MB, (
            f"Comprehensive test: memory grew by {final_stats.growth_mb:.2f}MB, "
            f"exceeds limit of {MAX_MEMORY_GROWTH_MB}MB"
        )
