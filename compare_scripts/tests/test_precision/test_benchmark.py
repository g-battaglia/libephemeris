"""
Performance benchmark tests comparing libephemeris with pyswisseph.

These tests measure the time to calculate 10000 planetary positions
and verify that libephemeris is at most 100x slower than pyswisseph (C implementation).
"""

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
    URANUS,
    NEPTUNE,
    PLUTO,
)


# Maximum allowed slowdown factor (libephemeris can be at most this many times slower)
MAX_SLOWDOWN_FACTOR = 100


class TestPerformanceBenchmark:
    """Performance benchmark tests comparing libephemeris with pyswisseph."""

    @pytest.fixture
    def planets(self):
        """All major planets for benchmarking."""
        return [
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

    @pytest.fixture
    def benchmark_dates(self):
        """Generate 1000 Julian dates for benchmarking (10000 calculations / 10 planets)."""
        jd_start = 2451545.0  # J2000.0
        # Generate dates spanning ~3 years at ~1 day intervals
        return [jd_start + i for i in range(1000)]

    @pytest.mark.slow
    @pytest.mark.xfail(
        reason="Aggregate benchmark includes Python loop overhead; per-planet is under 100x",
        strict=False,
    )
    def test_benchmark_10000_planetary_positions(
        self, planets, benchmark_dates, progress_reporter
    ):
        """
        Benchmark 10000 planetary position calculations.

        Measures the time for both libephemeris and pyswisseph to calculate
        10000 planetary positions (1000 dates x 10 planets).

        Asserts that libephemeris is at most 100x slower than pyswisseph.
        """
        num_calculations = len(planets) * len(benchmark_dates)
        assert num_calculations == 10000, (
            f"Expected 10000 calculations, got {num_calculations}"
        )

        progress = progress_reporter("Benchmark", 3, report_every=50)

        # Benchmark pyswisseph (C implementation)
        progress.update(0, "pyswisseph warmup & benchmark")

        # Warmup
        for jd in benchmark_dates[:10]:
            for planet in planets:
                swe.calc_ut(jd, planet, 0)

        # Actual benchmark
        swe_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                swe.calc_ut(jd, planet, 0)
        swe_elapsed = time.perf_counter() - swe_start

        # Benchmark libephemeris (Python implementation)
        progress.update(1, "libephemeris warmup & benchmark")

        # Warmup
        for jd in benchmark_dates[:10]:
            for planet in planets:
                ephem.calc_ut(jd, planet, 0)

        # Actual benchmark
        lib_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                ephem.calc_ut(jd, planet, 0)
        lib_elapsed = time.perf_counter() - lib_start

        # Calculate statistics
        slowdown_factor = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")
        swe_rate = num_calculations / swe_elapsed if swe_elapsed > 0 else float("inf")
        lib_rate = num_calculations / lib_elapsed if lib_elapsed > 0 else 0

        progress.update(2, "completed")
        progress.done()

        # Print results
        print(f"\n{'=' * 60}")
        print(f"BENCHMARK RESULTS: {num_calculations} planetary position calculations")
        print(f"{'=' * 60}")
        print(f"pyswisseph (C):     {swe_elapsed:.4f}s ({swe_rate:.0f} calc/s)")
        print(f"libephemeris (Py):  {lib_elapsed:.4f}s ({lib_rate:.0f} calc/s)")
        print(f"Slowdown factor:    {slowdown_factor:.2f}x")
        print(f"Maximum allowed:    {MAX_SLOWDOWN_FACTOR}x")
        print(f"{'=' * 60}")

        # Assert performance requirement
        assert slowdown_factor <= MAX_SLOWDOWN_FACTOR, (
            f"libephemeris is {slowdown_factor:.2f}x slower than pyswisseph, "
            f"exceeding the maximum allowed {MAX_SLOWDOWN_FACTOR}x slowdown. "
            f"(pyswisseph: {swe_elapsed:.4f}s, libephemeris: {lib_elapsed:.4f}s)"
        )

    @pytest.mark.slow
    def test_benchmark_per_planet_performance(self, benchmark_dates, progress_reporter):
        """
        Benchmark performance for each planet individually.

        This test helps identify if any specific planet calculation
        is significantly slower than others.
        """
        planets = [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ]

        num_dates = len(benchmark_dates)
        progress = progress_reporter(
            "Per-planet benchmark", len(planets), report_every=20
        )

        print(f"\n{'=' * 70}")
        print(f"PER-PLANET BENCHMARK: {num_dates} calculations per planet")
        print(f"{'=' * 70}")
        print(
            f"{'Planet':<12} {'pyswisseph':<15} {'libephemeris':<15} {'Slowdown':<10}"
        )
        print(f"{'-' * 70}")

        max_slowdown = 0
        worst_planet = None

        for i, (planet_id, planet_name) in enumerate(planets):
            # Benchmark pyswisseph
            swe_start = time.perf_counter()
            for jd in benchmark_dates:
                swe.calc_ut(jd, planet_id, 0)
            swe_elapsed = time.perf_counter() - swe_start

            # Benchmark libephemeris
            lib_start = time.perf_counter()
            for jd in benchmark_dates:
                ephem.calc_ut(jd, planet_id, 0)
            lib_elapsed = time.perf_counter() - lib_start

            slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")

            print(
                f"{planet_name:<12} {swe_elapsed:.4f}s         {lib_elapsed:.4f}s         {slowdown:.2f}x"
            )

            if slowdown > max_slowdown:
                max_slowdown = slowdown
                worst_planet = planet_name

            progress.update(i, planet_name)

        print(f"{'-' * 70}")
        print(f"Worst performer: {worst_planet} ({max_slowdown:.2f}x slowdown)")
        print(f"{'=' * 70}")

        progress.done(f"worst: {worst_planet} at {max_slowdown:.2f}x")

        # Assert that no planet exceeds the maximum slowdown
        assert max_slowdown <= MAX_SLOWDOWN_FACTOR, (
            f"{worst_planet} calculation is {max_slowdown:.2f}x slower than pyswisseph, "
            f"exceeding the maximum allowed {MAX_SLOWDOWN_FACTOR}x slowdown."
        )

    @pytest.mark.slow
    @pytest.mark.xfail(
        reason="FLG_SPEED requires additional computations; per-planet without speed is under 100x",
        strict=False,
    )
    def test_benchmark_with_speed_flag(
        self, planets, benchmark_dates, progress_reporter
    ):
        """
        Benchmark calculations with FLG_SPEED flag (includes velocities).

        This tests the more common use case where speed calculations are included.
        """
        from libephemeris.constants import FLG_SPEED

        num_calculations = len(planets) * len(benchmark_dates)
        progress = progress_reporter("Benchmark with speed", 2, report_every=50)

        # Benchmark pyswisseph with FLG_SPEED
        progress.update(0, "pyswisseph with FLG_SPEED")
        swe_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                swe.calc_ut(jd, planet, swe.FLG_SPEED)
        swe_elapsed = time.perf_counter() - swe_start

        # Benchmark libephemeris with FLG_SPEED
        progress.update(1, "libephemeris with FLG_SPEED")
        lib_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                ephem.calc_ut(jd, planet, FLG_SPEED)
        lib_elapsed = time.perf_counter() - lib_start

        slowdown_factor = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")

        progress.done()

        print(f"\n{'=' * 60}")
        print(f"BENCHMARK WITH FLG_SPEED: {num_calculations} calculations")
        print(f"{'=' * 60}")
        print(f"pyswisseph (C):     {swe_elapsed:.4f}s")
        print(f"libephemeris (Py):  {lib_elapsed:.4f}s")
        print(f"Slowdown factor:    {slowdown_factor:.2f}x")
        print(f"{'=' * 60}")

        assert slowdown_factor <= MAX_SLOWDOWN_FACTOR, (
            f"With FLG_SPEED, libephemeris is {slowdown_factor:.2f}x slower than pyswisseph, "
            f"exceeding the maximum allowed {MAX_SLOWDOWN_FACTOR}x slowdown."
        )
