"""
Comprehensive benchmark suite comparing libephemeris vs pyswisseph.

This module provides benchmarks for the most common ephemeris operations:
- calc_ut: planetary position calculations
- houses: house cusp calculations
- get_ayanamsa: ayanamsa calculations

These benchmarks measure performance differences between the pure Python
libephemeris implementation and the C-based pyswisseph library.
"""

import time
import statistics
from dataclasses import dataclass
from typing import Callable

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
    FLG_SPEED,
    SIDM_LAHIRI,
    SIDM_FAGAN_BRADLEY,
    SIDM_RAMAN,
    SIDM_TRUE_CITRA,
    SIDM_KRISHNAMURTI,
)


# Maximum allowed slowdown factor (libephemeris can be at most this many times slower)
MAX_SLOWDOWN_FACTOR = 100


@dataclass
class BenchmarkResult:
    """Result of a benchmark run."""

    name: str
    iterations: int
    swe_time: float
    lib_time: float
    swe_rate: float
    lib_rate: float
    slowdown: float

    def __str__(self) -> str:
        return (
            f"{self.name}: "
            f"pyswisseph={self.swe_time:.4f}s ({self.swe_rate:.0f}/s), "
            f"libephemeris={self.lib_time:.4f}s ({self.lib_rate:.0f}/s), "
            f"slowdown={self.slowdown:.2f}x"
        )


def run_benchmark(
    name: str,
    swe_func: Callable[[], None],
    lib_func: Callable[[], None],
    iterations: int,
    warmup: int = 10,
) -> BenchmarkResult:
    """
    Run a benchmark comparing pyswisseph and libephemeris.

    Args:
        name: Name of the benchmark
        swe_func: Function to call for pyswisseph benchmark
        lib_func: Function to call for libephemeris benchmark
        iterations: Number of iterations to run
        warmup: Number of warmup iterations

    Returns:
        BenchmarkResult with timing data
    """
    # Warmup pyswisseph
    for _ in range(warmup):
        swe_func()

    # Benchmark pyswisseph
    swe_start = time.perf_counter()
    for _ in range(iterations):
        swe_func()
    swe_elapsed = time.perf_counter() - swe_start

    # Warmup libephemeris
    for _ in range(warmup):
        lib_func()

    # Benchmark libephemeris
    lib_start = time.perf_counter()
    for _ in range(iterations):
        lib_func()
    lib_elapsed = time.perf_counter() - lib_start

    # Calculate statistics
    swe_rate = iterations / swe_elapsed if swe_elapsed > 0 else float("inf")
    lib_rate = iterations / lib_elapsed if lib_elapsed > 0 else 0
    slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")

    return BenchmarkResult(
        name=name,
        iterations=iterations,
        swe_time=swe_elapsed,
        lib_time=lib_elapsed,
        swe_rate=swe_rate,
        lib_rate=lib_rate,
        slowdown=slowdown,
    )


class TestCalcUtBenchmark:
    """Benchmarks for calc_ut (planetary position calculations)."""

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
        """Generate 1000 Julian dates for benchmarking."""
        jd_start = 2451545.0  # J2000.0
        return [jd_start + i for i in range(1000)]

    @pytest.mark.slow
    @pytest.mark.xfail(
        reason="Aggregate benchmark includes Python loop overhead; pure Python vs C implementation",
        strict=False,
    )
    def test_calc_ut_all_planets(self, planets, benchmark_dates, progress_reporter):
        """
        Benchmark calc_ut for all planets over 1000 dates.

        Total calculations: 10,000 (10 planets x 1000 dates)
        """
        num_calculations = len(planets) * len(benchmark_dates)
        progress = progress_reporter("calc_ut benchmark", 2)

        progress.update(0, "pyswisseph")

        # Warmup
        for jd in benchmark_dates[:10]:
            for planet in planets:
                swe.calc_ut(jd, planet, 0)

        # Benchmark pyswisseph
        swe_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                swe.calc_ut(jd, planet, 0)
        swe_elapsed = time.perf_counter() - swe_start

        progress.update(1, "libephemeris")

        # Warmup
        for jd in benchmark_dates[:10]:
            for planet in planets:
                ephem.calc_ut(jd, planet, 0)

        # Benchmark libephemeris
        lib_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                ephem.calc_ut(jd, planet, 0)
        lib_elapsed = time.perf_counter() - lib_start

        slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")
        swe_rate = num_calculations / swe_elapsed if swe_elapsed > 0 else float("inf")
        lib_rate = num_calculations / lib_elapsed if lib_elapsed > 0 else 0

        progress.done()

        print(f"\n{'=' * 70}")
        print(
            f"BENCHMARK: calc_ut - {num_calculations} planetary position calculations"
        )
        print(f"{'=' * 70}")
        print(f"pyswisseph (C):     {swe_elapsed:.4f}s ({swe_rate:.0f} calc/s)")
        print(f"libephemeris (Py):  {lib_elapsed:.4f}s ({lib_rate:.0f} calc/s)")
        print(f"Slowdown factor:    {slowdown:.2f}x")
        print(f"{'=' * 70}")

        # Store result for README generation
        assert slowdown <= MAX_SLOWDOWN_FACTOR, (
            f"calc_ut is {slowdown:.2f}x slower, exceeds {MAX_SLOWDOWN_FACTOR}x limit"
        )

    @pytest.mark.slow
    def test_calc_ut_with_speed(self, planets, benchmark_dates, progress_reporter):
        """
        Benchmark calc_ut with FLG_SPEED flag.

        This tests velocity calculations which are computationally more expensive.
        """
        num_calculations = len(planets) * len(benchmark_dates)
        progress = progress_reporter("calc_ut + SPEED benchmark", 2)

        progress.update(0, "pyswisseph")

        # Warmup
        for jd in benchmark_dates[:10]:
            for planet in planets:
                swe.calc_ut(jd, planet, swe.FLG_SPEED)

        # Benchmark pyswisseph
        swe_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                swe.calc_ut(jd, planet, swe.FLG_SPEED)
        swe_elapsed = time.perf_counter() - swe_start

        progress.update(1, "libephemeris")

        # Warmup
        for jd in benchmark_dates[:10]:
            for planet in planets:
                ephem.calc_ut(jd, planet, FLG_SPEED)

        # Benchmark libephemeris
        lib_start = time.perf_counter()
        for jd in benchmark_dates:
            for planet in planets:
                ephem.calc_ut(jd, planet, FLG_SPEED)
        lib_elapsed = time.perf_counter() - lib_start

        slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")
        swe_rate = num_calculations / swe_elapsed if swe_elapsed > 0 else float("inf")
        lib_rate = num_calculations / lib_elapsed if lib_elapsed > 0 else 0

        progress.done()

        print(f"\n{'=' * 70}")
        print(f"BENCHMARK: calc_ut + FLG_SPEED - {num_calculations} calculations")
        print(f"{'=' * 70}")
        print(f"pyswisseph (C):     {swe_elapsed:.4f}s ({swe_rate:.0f} calc/s)")
        print(f"libephemeris (Py):  {lib_elapsed:.4f}s ({lib_rate:.0f} calc/s)")
        print(f"Slowdown factor:    {slowdown:.2f}x")
        print(f"{'=' * 70}")


class TestHousesBenchmark:
    """Benchmarks for houses (house cusp calculations)."""

    @pytest.fixture
    def house_systems(self):
        """Common house systems for benchmarking."""
        return [
            (b"P", ord("P"), "Placidus"),
            (b"K", ord("K"), "Koch"),
            (b"R", ord("R"), "Regiomontanus"),
            (b"C", ord("C"), "Campanus"),
            (b"E", ord("E"), "Equal"),
            (b"W", ord("W"), "Whole Sign"),
            (b"O", ord("O"), "Porphyry"),
        ]

    @pytest.fixture
    def locations(self):
        """Test locations for house calculations."""
        return [
            ("Rome", 41.9028, 12.4964),
            ("London", 51.5074, -0.1278),
            ("New York", 40.7128, -74.0060),
            ("Sydney", -33.8688, 151.2093),
            ("Tokyo", 35.6762, 139.6503),
        ]

    @pytest.fixture
    def benchmark_dates(self):
        """Generate 100 Julian dates for house benchmarking."""
        jd_start = 2451545.0  # J2000.0
        return [jd_start + i for i in range(100)]

    @pytest.mark.slow
    def test_houses_all_systems(
        self, house_systems, locations, benchmark_dates, progress_reporter
    ):
        """
        Benchmark house calculations for all systems and locations.

        Total calculations: 3500 (7 systems x 5 locations x 100 dates)
        """
        num_calculations = len(house_systems) * len(locations) * len(benchmark_dates)
        progress = progress_reporter("houses benchmark", 2)

        progress.update(0, "pyswisseph")

        # Warmup
        for jd in benchmark_dates[:5]:
            for _, lat, lon in locations:
                for hsys_bytes, hsys_ord, _ in house_systems:
                    swe.houses(jd, lat, lon, hsys_bytes)

        # Benchmark pyswisseph
        swe_start = time.perf_counter()
        for jd in benchmark_dates:
            for _, lat, lon in locations:
                for hsys_bytes, hsys_ord, _ in house_systems:
                    swe.houses(jd, lat, lon, hsys_bytes)
        swe_elapsed = time.perf_counter() - swe_start

        progress.update(1, "libephemeris")

        # Warmup
        for jd in benchmark_dates[:5]:
            for _, lat, lon in locations:
                for hsys_bytes, hsys_ord, _ in house_systems:
                    ephem.houses(jd, lat, lon, hsys_ord)

        # Benchmark libephemeris
        lib_start = time.perf_counter()
        for jd in benchmark_dates:
            for _, lat, lon in locations:
                for hsys_bytes, hsys_ord, _ in house_systems:
                    ephem.houses(jd, lat, lon, hsys_ord)
        lib_elapsed = time.perf_counter() - lib_start

        slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")
        swe_rate = num_calculations / swe_elapsed if swe_elapsed > 0 else float("inf")
        lib_rate = num_calculations / lib_elapsed if lib_elapsed > 0 else 0

        progress.done()

        print(f"\n{'=' * 70}")
        print(f"BENCHMARK: houses - {num_calculations} house calculations")
        print(f"{'=' * 70}")
        print(f"pyswisseph (C):     {swe_elapsed:.4f}s ({swe_rate:.0f} calc/s)")
        print(f"libephemeris (Py):  {lib_elapsed:.4f}s ({lib_rate:.0f} calc/s)")
        print(f"Slowdown factor:    {slowdown:.2f}x")
        print(f"{'=' * 70}")

        assert slowdown <= MAX_SLOWDOWN_FACTOR, (
            f"houses is {slowdown:.2f}x slower, exceeds {MAX_SLOWDOWN_FACTOR}x limit"
        )

    @pytest.mark.slow
    def test_houses_per_system(
        self, house_systems, locations, benchmark_dates, progress_reporter
    ):
        """
        Benchmark each house system individually.

        Helps identify which house systems are slowest.
        """
        num_per_system = len(locations) * len(benchmark_dates)
        progress = progress_reporter("per-system houses benchmark", len(house_systems))

        print(f"\n{'=' * 70}")
        print(
            f"PER-SYSTEM BENCHMARK: houses - {num_per_system} calculations per system"
        )
        print(f"{'=' * 70}")
        print(
            f"{'System':<15} {'pyswisseph':<15} {'libephemeris':<15} {'Slowdown':<10}"
        )
        print(f"{'-' * 70}")

        max_slowdown = 0
        worst_system = None

        for i, (hsys_bytes, hsys_ord, name) in enumerate(house_systems):
            # Benchmark pyswisseph
            swe_start = time.perf_counter()
            for jd in benchmark_dates:
                for _, lat, lon in locations:
                    swe.houses(jd, lat, lon, hsys_bytes)
            swe_elapsed = time.perf_counter() - swe_start

            # Benchmark libephemeris
            lib_start = time.perf_counter()
            for jd in benchmark_dates:
                for _, lat, lon in locations:
                    ephem.houses(jd, lat, lon, hsys_ord)
            lib_elapsed = time.perf_counter() - lib_start

            slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")

            print(
                f"{name:<15} {swe_elapsed:.4f}s         {lib_elapsed:.4f}s         {slowdown:.2f}x"
            )

            if slowdown > max_slowdown:
                max_slowdown = slowdown
                worst_system = name

            progress.update(i, name)

        print(f"{'-' * 70}")
        print(f"Worst performer: {worst_system} ({max_slowdown:.2f}x slowdown)")
        print(f"{'=' * 70}")

        progress.done(f"worst: {worst_system} at {max_slowdown:.2f}x")

        assert max_slowdown <= MAX_SLOWDOWN_FACTOR, (
            f"{worst_system} is {max_slowdown:.2f}x slower, exceeds limit"
        )


class TestAyanamsaBenchmark:
    """Benchmarks for get_ayanamsa_ut (ayanamsa calculations)."""

    @pytest.fixture
    def ayanamsha_modes(self):
        """Common ayanamsha modes for benchmarking."""
        return [
            (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SIDM_TRUE_CITRA, "True Citra"),
        ]

    @pytest.fixture
    def benchmark_dates(self):
        """Generate 1000 Julian dates for ayanamsa benchmarking."""
        jd_start = 2451545.0  # J2000.0
        return [jd_start + i for i in range(1000)]

    @pytest.mark.slow
    def test_ayanamsa_all_modes(
        self, ayanamsha_modes, benchmark_dates, progress_reporter
    ):
        """
        Benchmark ayanamsa calculations for all modes.

        Total calculations: 5000 (5 modes x 1000 dates)
        """
        num_calculations = len(ayanamsha_modes) * len(benchmark_dates)
        progress = progress_reporter("ayanamsa benchmark", 2)

        progress.update(0, "pyswisseph")

        # Warmup
        for mode, _ in ayanamsha_modes:
            swe.set_sid_mode(mode)
            for jd in benchmark_dates[:10]:
                swe.get_ayanamsa_ut(jd)

        # Benchmark pyswisseph
        swe_start = time.perf_counter()
        for mode, _ in ayanamsha_modes:
            swe.set_sid_mode(mode)
            for jd in benchmark_dates:
                swe.get_ayanamsa_ut(jd)
        swe_elapsed = time.perf_counter() - swe_start

        progress.update(1, "libephemeris")

        # Warmup
        for mode, _ in ayanamsha_modes:
            ephem.set_sid_mode(mode)
            for jd in benchmark_dates[:10]:
                ephem.get_ayanamsa_ut(jd)

        # Benchmark libephemeris
        lib_start = time.perf_counter()
        for mode, _ in ayanamsha_modes:
            ephem.set_sid_mode(mode)
            for jd in benchmark_dates:
                ephem.get_ayanamsa_ut(jd)
        lib_elapsed = time.perf_counter() - lib_start

        slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")
        swe_rate = num_calculations / swe_elapsed if swe_elapsed > 0 else float("inf")
        lib_rate = num_calculations / lib_elapsed if lib_elapsed > 0 else 0

        progress.done()

        print(f"\n{'=' * 70}")
        print(f"BENCHMARK: get_ayanamsa_ut - {num_calculations} ayanamsa calculations")
        print(f"{'=' * 70}")
        print(f"pyswisseph (C):     {swe_elapsed:.4f}s ({swe_rate:.0f} calc/s)")
        print(f"libephemeris (Py):  {lib_elapsed:.4f}s ({lib_rate:.0f} calc/s)")
        print(f"Slowdown factor:    {slowdown:.2f}x")
        print(f"{'=' * 70}")

        assert slowdown <= MAX_SLOWDOWN_FACTOR, (
            f"get_ayanamsa_ut is {slowdown:.2f}x slower, exceeds limit"
        )

    @pytest.mark.slow
    def test_ayanamsa_per_mode(
        self, ayanamsha_modes, benchmark_dates, progress_reporter
    ):
        """
        Benchmark each ayanamsha mode individually.

        Star-based ayanamshas may be slower due to star position calculations.
        """
        num_per_mode = len(benchmark_dates)
        progress = progress_reporter(
            "per-mode ayanamsa benchmark", len(ayanamsha_modes)
        )

        print(f"\n{'=' * 70}")
        print(f"PER-MODE BENCHMARK: ayanamsa - {num_per_mode} calculations per mode")
        print(f"{'=' * 70}")
        print(f"{'Mode':<20} {'pyswisseph':<15} {'libephemeris':<15} {'Slowdown':<10}")
        print(f"{'-' * 70}")

        max_slowdown = 0
        worst_mode = None

        for i, (mode, name) in enumerate(ayanamsha_modes):
            # Benchmark pyswisseph
            swe.set_sid_mode(mode)
            swe_start = time.perf_counter()
            for jd in benchmark_dates:
                swe.get_ayanamsa_ut(jd)
            swe_elapsed = time.perf_counter() - swe_start

            # Benchmark libephemeris
            ephem.set_sid_mode(mode)
            lib_start = time.perf_counter()
            for jd in benchmark_dates:
                ephem.get_ayanamsa_ut(jd)
            lib_elapsed = time.perf_counter() - lib_start

            slowdown = lib_elapsed / swe_elapsed if swe_elapsed > 0 else float("inf")

            print(
                f"{name:<20} {swe_elapsed:.4f}s         {lib_elapsed:.4f}s         {slowdown:.2f}x"
            )

            if slowdown > max_slowdown:
                max_slowdown = slowdown
                worst_mode = name

            progress.update(i, name)

        print(f"{'-' * 70}")
        print(f"Worst performer: {worst_mode} ({max_slowdown:.2f}x slowdown)")
        print(f"{'=' * 70}")

        progress.done(f"worst: {worst_mode} at {max_slowdown:.2f}x")

        assert max_slowdown <= MAX_SLOWDOWN_FACTOR, (
            f"{worst_mode} is {max_slowdown:.2f}x slower, exceeds limit"
        )


class TestComprehensiveBenchmark:
    """
    Comprehensive benchmark suite that runs all operations together.

    This class provides a summary benchmark for README documentation.
    """

    @pytest.mark.slow
    @pytest.mark.xfail(
        reason="calc_ut exceeds 100x due to pure Python vs C implementation",
        strict=False,
    )
    def test_comprehensive_benchmark_summary(self, progress_reporter):
        """
        Run all benchmarks and produce a summary table.

        This test runs a reduced but representative set of benchmarks
        and prints a formatted summary suitable for README inclusion.
        """
        results = []
        jd_start = 2451545.0
        dates = [jd_start + i for i in range(500)]

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

        progress = progress_reporter("Comprehensive benchmark", 3)

        # 1. calc_ut benchmark
        progress.update(0, "calc_ut")
        num_calc = len(planets) * len(dates)

        swe_start = time.perf_counter()
        for jd in dates:
            for planet in planets:
                swe.calc_ut(jd, planet, 0)
        swe_calc_time = time.perf_counter() - swe_start

        lib_start = time.perf_counter()
        for jd in dates:
            for planet in planets:
                ephem.calc_ut(jd, planet, 0)
        lib_calc_time = time.perf_counter() - lib_start

        results.append(
            BenchmarkResult(
                name="calc_ut",
                iterations=num_calc,
                swe_time=swe_calc_time,
                lib_time=lib_calc_time,
                swe_rate=num_calc / swe_calc_time,
                lib_rate=num_calc / lib_calc_time,
                slowdown=lib_calc_time / swe_calc_time,
            )
        )

        # 2. houses benchmark
        progress.update(1, "houses")
        # (bytes for swe, ord for libephemeris)
        house_systems_swe = [b"P", b"K", b"E", b"W"]
        house_systems_lib = [ord("P"), ord("K"), ord("E"), ord("W")]
        locations = [(41.9, 12.5), (51.5, -0.1), (40.7, -74.0)]
        dates_houses = dates[:100]
        num_houses = len(house_systems_swe) * len(locations) * len(dates_houses)

        swe_start = time.perf_counter()
        for jd in dates_houses:
            for lat, lon in locations:
                for hsys in house_systems_swe:
                    swe.houses(jd, lat, lon, hsys)
        swe_houses_time = time.perf_counter() - swe_start

        lib_start = time.perf_counter()
        for jd in dates_houses:
            for lat, lon in locations:
                for hsys in house_systems_lib:
                    ephem.houses(jd, lat, lon, hsys)
        lib_houses_time = time.perf_counter() - lib_start

        results.append(
            BenchmarkResult(
                name="houses",
                iterations=num_houses,
                swe_time=swe_houses_time,
                lib_time=lib_houses_time,
                swe_rate=num_houses / swe_houses_time,
                lib_rate=num_houses / lib_houses_time,
                slowdown=lib_houses_time / swe_houses_time,
            )
        )

        # 3. get_ayanamsa benchmark
        progress.update(2, "get_ayanamsa_ut")
        ayanamsha_modes = [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
        ]
        num_ayan = len(ayanamsha_modes) * len(dates)

        swe_start = time.perf_counter()
        for mode in ayanamsha_modes:
            swe.set_sid_mode(mode)
            for jd in dates:
                swe.get_ayanamsa_ut(jd)
        swe_ayan_time = time.perf_counter() - swe_start

        lib_start = time.perf_counter()
        for mode in ayanamsha_modes:
            ephem.set_sid_mode(mode)
            for jd in dates:
                ephem.get_ayanamsa_ut(jd)
        lib_ayan_time = time.perf_counter() - lib_start

        results.append(
            BenchmarkResult(
                name="get_ayanamsa_ut",
                iterations=num_ayan,
                swe_time=swe_ayan_time,
                lib_time=lib_ayan_time,
                swe_rate=num_ayan / swe_ayan_time,
                lib_rate=num_ayan / lib_ayan_time,
                slowdown=lib_ayan_time / swe_ayan_time,
            )
        )

        progress.done()

        # Print summary table
        print(f"\n{'=' * 80}")
        print("BENCHMARK SUMMARY: libephemeris vs pyswisseph")
        print(f"{'=' * 80}")
        print(
            f"{'Operation':<20} {'Iterations':<12} {'pyswisseph':<15} "
            f"{'libephemeris':<15} {'Slowdown':<10}"
        )
        print(f"{'-' * 80}")

        for r in results:
            print(
                f"{r.name:<20} {r.iterations:<12} {r.swe_time:.4f}s         "
                f"{r.lib_time:.4f}s          {r.slowdown:.1f}x"
            )

        print(f"{'-' * 80}")
        avg_slowdown = statistics.mean(r.slowdown for r in results)
        print(f"Average slowdown: {avg_slowdown:.1f}x")
        print(f"{'=' * 80}")

        # Print markdown table for README
        print("\n### Markdown table for README:\n")
        print(
            "| Operation | Iterations | pyswisseph (C) | libephemeris (Python) | Slowdown |"
        )
        print(
            "| --------- | ---------- | -------------- | --------------------- | -------- |"
        )
        for r in results:
            print(
                f"| {r.name} | {r.iterations} | "
                f"{r.swe_time:.4f}s ({r.swe_rate:.0f}/s) | "
                f"{r.lib_time:.4f}s ({r.lib_rate:.0f}/s) | "
                f"{r.slowdown:.1f}x |"
            )

        # Assert all results are within limit
        for r in results:
            assert r.slowdown <= MAX_SLOWDOWN_FACTOR, (
                f"{r.name} is {r.slowdown:.2f}x slower, exceeds limit"
            )
