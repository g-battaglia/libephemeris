"""Performance benchmark tests for regression detection.

Establishes baseline timings for key operations and alerts on significant
regressions. Uses manual timing with perf_counter.

Validation Plan v2, Section 7.2.

Targets (module-level API):
    - calc_ut (Skyfield): < 5000 µs (first-call overhead)
    - calc_ut (LEB): < 500 µs
    - houses: < 500 µs
    - LEB speedup: >= 8x over Skyfield for JPL-backed planets

Note: The module-level API has overhead from the context-swap lock and
global state save/restore. The EphemerisContext API avoids this overhead.
Absolute timings vary by platform; the key metric is the LEB/Skyfield
speedup ratio.
"""

from __future__ import annotations

import os
import time
import warnings

import pytest

import libephemeris as swe

warnings.filterwarnings("ignore")

pytestmark = pytest.mark.slow

J2000 = 2451545.0

# Number of iterations for timing
WARMUP = 50
ITERATIONS = 500


def measure_us(func, warmup: int = WARMUP, iters: int = ITERATIONS) -> float:
    """Measure mean execution time in microseconds."""
    for _ in range(warmup):
        func()

    start = time.perf_counter()
    for _ in range(iters):
        func()
    elapsed = time.perf_counter() - start

    return (elapsed / iters) * 1e6


class TestPerformanceBaselines:
    """§7.2 Performance benchmarks for key operations."""

    def test_calc_ut_performance(self) -> None:
        """calc_ut should complete in reasonable time."""

        def calc():
            swe.calc_ut(J2000, swe.SUN, swe.FLG_SPEED)

        us = measure_us(calc)
        print(f"\n  calc_ut (Sun, SPEED): {us:.1f} us")
        assert us < 5000, f"calc_ut took {us:.1f} us (hard limit 5000 us)"

    def test_calc_ut_moon_performance(self) -> None:
        """calc_ut Moon performance."""

        def calc():
            swe.calc_ut(J2000, swe.MOON, swe.FLG_SPEED)

        us = measure_us(calc)
        print(f"\n  calc_ut (Moon, SPEED): {us:.1f} us")
        assert us < 5000, f"calc_ut Moon took {us:.1f} us (hard limit 5000 us)"

    def test_houses_performance(self) -> None:
        """houses() should be < 500 us per call."""

        def calc():
            swe.houses(J2000, 41.9, 12.5, ord("P"))

        us = measure_us(calc)
        print(f"\n  houses (Placidus): {us:.1f} us")
        assert us < 2000, f"houses took {us:.1f} us (hard limit 2000 us)"

    def test_houses_equal_performance(self) -> None:
        """Equal house system performance."""

        def calc():
            swe.houses(J2000, 41.9, 12.5, ord("E"))

        us = measure_us(calc)
        print(f"\n  houses (Equal): {us:.1f} us")
        assert us < 2000, f"houses Equal took {us:.1f} us (hard limit 2000 us)"

    def test_sidtime_performance(self) -> None:
        """sidtime() should be very fast."""

        def calc():
            swe.sidtime(J2000)

        us = measure_us(calc, iters=1000)
        print(f"\n  sidtime: {us:.1f} us")
        assert us < 200, f"sidtime took {us:.1f} us (hard limit 200 us)"

    def test_deltat_performance(self) -> None:
        """deltat() should be very fast."""

        def calc():
            swe.deltat(J2000)

        us = measure_us(calc, iters=1000)
        print(f"\n  deltat: {us:.1f} us")
        assert us < 200, f"deltat took {us:.1f} us (hard limit 200 us)"


class TestLEBSpeedupRatio:
    """§7.2 LEB vs Skyfield speedup ratio — the key performance metric."""

    def _find_leb_path(self) -> str:
        """Find an available LEB file or skip."""
        for candidate in [
            "data/leb/ephemeris_medium.leb",
            "data/leb/ephemeris_base.leb",
        ]:
            if os.path.exists(candidate):
                return candidate
        pytest.skip("No LEB file available")
        return ""  # unreachable

    def test_leb_speedup_sun(self) -> None:
        """LEB should provide significant speedup for Sun."""
        leb_path = self._find_leb_path()

        # Measure Skyfield mode
        swe.set_calc_mode("skyfield")

        def calc_sky():
            swe.calc_ut(J2000, swe.SUN, swe.FLG_SPEED)

        us_sky = measure_us(calc_sky)

        # Measure LEB mode
        swe.set_calc_mode("leb")
        swe.set_leb_file(leb_path)

        def calc_leb():
            swe.calc_ut(J2000, swe.SUN, swe.FLG_SPEED)

        us_leb = measure_us(calc_leb)

        # Restore default
        swe.set_calc_mode(None)

        ratio = us_sky / us_leb if us_leb > 0 else float("inf")
        print(
            f"\n  LEB speedup (Sun): {ratio:.1f}x "
            f"(Skyfield={us_sky:.1f}us, LEB={us_leb:.1f}us)"
        )
        assert ratio >= 4.0, (
            f"LEB speedup only {ratio:.1f}x for Sun (target >= 8x, hard limit >= 4x)"
        )

    def test_leb_speedup_mars(self) -> None:
        """LEB should provide >= 8x speedup for Mars."""
        leb_path = self._find_leb_path()

        # Measure Skyfield mode
        swe.set_calc_mode("skyfield")

        def calc_sky():
            swe.calc_ut(J2000, swe.MARS, swe.FLG_SPEED)

        us_sky = measure_us(calc_sky)

        # Measure LEB mode
        swe.set_calc_mode("leb")
        swe.set_leb_file(leb_path)

        def calc_leb():
            swe.calc_ut(J2000, swe.MARS, swe.FLG_SPEED)

        us_leb = measure_us(calc_leb)

        # Restore default
        swe.set_calc_mode(None)

        ratio = us_sky / us_leb if us_leb > 0 else float("inf")
        print(
            f"\n  LEB speedup (Mars): {ratio:.1f}x "
            f"(Skyfield={us_sky:.1f}us, LEB={us_leb:.1f}us)"
        )
        assert ratio >= 4.0, (
            f"LEB speedup only {ratio:.1f}x for Mars (target >= 8x, hard limit >= 4x)"
        )

    def test_leb_speedup_moon(self) -> None:
        """LEB speedup for Moon (different code path)."""
        leb_path = self._find_leb_path()

        # Measure Skyfield mode
        swe.set_calc_mode("skyfield")

        def calc_sky():
            swe.calc_ut(J2000, swe.MOON, swe.FLG_SPEED)

        us_sky = measure_us(calc_sky)

        # Measure LEB mode
        swe.set_calc_mode("leb")
        swe.set_leb_file(leb_path)

        def calc_leb():
            swe.calc_ut(J2000, swe.MOON, swe.FLG_SPEED)

        us_leb = measure_us(calc_leb)

        # Restore default
        swe.set_calc_mode(None)

        ratio = us_sky / us_leb if us_leb > 0 else float("inf")
        print(
            f"\n  LEB speedup (Moon): {ratio:.1f}x "
            f"(Skyfield={us_sky:.1f}us, LEB={us_leb:.1f}us)"
        )
        # Moon may have lower speedup due to different calculation path
        assert ratio >= 2.0, (
            f"LEB speedup only {ratio:.1f}x for Moon (hard limit >= 2x)"
        )


class TestPerformanceSummary:
    """Print a comprehensive performance summary."""

    def test_print_summary(self) -> None:
        """Collect and print a performance summary table."""
        leb_path = None
        for candidate in [
            "data/leb/ephemeris_medium.leb",
            "data/leb/ephemeris_base.leb",
        ]:
            if os.path.exists(candidate):
                leb_path = candidate
                break

        results = []

        # Skyfield calc_ut (use explicit mode for accurate comparison)
        swe.set_calc_mode("skyfield")
        for body, name in [
            (swe.SUN, "Sun"),
            (swe.MOON, "Moon"),
            (swe.MARS, "Mars"),
            (swe.JUPITER, "Jupiter"),
        ]:

            def calc(b=body):
                swe.calc_ut(J2000, b, swe.FLG_SPEED)

            us = measure_us(calc)
            results.append(("calc_ut Skyfield", name, us))

        # LEB calc_ut
        if leb_path:
            swe.set_calc_mode("leb")
            swe.set_leb_file(leb_path)
            for body, name in [
                (swe.SUN, "Sun"),
                (swe.MOON, "Moon"),
                (swe.MARS, "Mars"),
                (swe.JUPITER, "Jupiter"),
            ]:

                def calc(b=body):
                    swe.calc_ut(J2000, b, swe.FLG_SPEED)

                us = measure_us(calc)
                results.append(("calc_ut LEB", name, us))

        # Restore default mode
        swe.set_calc_mode(None)

        # Houses
        for hsys, name in [
            (ord("P"), "Placidus"),
            (ord("E"), "Equal"),
            (ord("W"), "WholeSgn"),
        ]:

            def calc(h=hsys):
                swe.houses(J2000, 41.9, 12.5, h)

            us = measure_us(calc)
            results.append(("houses", name, us))

        # Utility functions
        def st():
            swe.sidtime(J2000)

        results.append(("sidtime", "-", measure_us(st, iters=1000)))

        def dt():
            swe.deltat(J2000)

        results.append(("deltat", "-", measure_us(dt, iters=1000)))

        # Print summary
        print("\n\n  Performance Summary")
        print("  " + "-" * 50)
        print(f"  {'Operation':<25} {'Body/Sys':<12} {'us':>8}")
        print("  " + "-" * 50)
        for op, detail, us in results:
            print(f"  {op:<25} {detail:<12} {us:>8.1f}")
        print("  " + "-" * 50)

        # Print speedup ratios if LEB data available
        if leb_path:
            print("\n  Speedup Ratios (Skyfield / LEB)")
            print("  " + "-" * 35)
            sky_results = {r[1]: r[2] for r in results if r[0] == "calc_ut Skyfield"}
            leb_results = {r[1]: r[2] for r in results if r[0] == "calc_ut LEB"}
            for name in sky_results:
                if name in leb_results and leb_results[name] > 0:
                    ratio = sky_results[name] / leb_results[name]
                    print(f"  {name:<12} {ratio:>6.1f}x")
            print("  " + "-" * 35)

        assert True  # Summary always passes
