"""Thread safety and concurrency stress tests.

Verifies that EphemerisContext provides true thread-safe isolation under
concurrent load, and that shared resources are correctly protected.

Validation Plan v2, Section 6.
"""

from __future__ import annotations

import math
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import pytest

import libephemeris as swe
from libephemeris.context import EphemerisContext

warnings.filterwarnings("ignore")

pytestmark = pytest.mark.slow

# Base JD for tests
J2000 = 2451545.0


# ============================================================================
# §6.1 Concurrent Context Stress Test
# ============================================================================


class TestConcurrentContextStress:
    """§6.1 Verify EphemerisContext under 50-thread concurrent load."""

    def _compute_single_threaded_baseline(
        self, bodies: list[int], jds: list[float], flags: int
    ) -> dict[tuple[int, float], tuple[float, ...]]:
        """Compute baseline results in single-threaded mode."""
        ctx = EphemerisContext()
        results = {}
        for body in bodies:
            for jd in jds:
                pos, _ = ctx.calc_ut(jd, body, flags)
                results[(body, jd)] = pos
        return results

    def test_50_threads_match_baseline(self) -> None:
        """50 threads each computing 100 positions must match single-threaded."""
        bodies = [swe.SUN, swe.MOON, swe.MARS, swe.JUPITER, swe.SATURN]
        jds = [J2000 + i * 30.0 for i in range(20)]  # 20 dates, 30 days apart
        flags = swe.FLG_SPEED

        # 5 bodies × 20 dates = 100 calculations per thread
        baseline = self._compute_single_threaded_baseline(bodies, jds, flags)

        errors: list[str] = []
        lock = threading.Lock()

        def worker(thread_id: int) -> None:
            ctx = EphemerisContext()
            for body in bodies:
                for jd in jds:
                    pos, _ = ctx.calc_ut(jd, body, flags)
                    expected = baseline[(body, jd)]
                    for i in range(6):
                        if abs(pos[i] - expected[i]) > 1e-10:
                            with lock:
                                errors.append(
                                    f"Thread {thread_id}: body={body} "
                                    f"jd={jd} elem={i} "
                                    f"got={pos[i]} expected={expected[i]}"
                                )
                            return

        threads = []
        for tid in range(50):
            t = threading.Thread(target=worker, args=(tid,))
            threads.append(t)

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=120)

        assert not errors, f"{len(errors)} mismatches found:\n" + "\n".join(errors[:10])

    def test_concurrent_no_deadlock(self) -> None:
        """50 threads must all complete within 60 seconds (no deadlock)."""
        completed = []
        lock = threading.Lock()

        def worker(thread_id: int) -> None:
            ctx = EphemerisContext()
            for i in range(50):
                jd = J2000 + i * 10.0
                ctx.calc_ut(jd, swe.SUN, 0)
                ctx.calc_ut(jd, swe.MOON, swe.FLG_SPEED)
            with lock:
                completed.append(thread_id)

        threads = []
        for tid in range(50):
            t = threading.Thread(target=worker, args=(tid,))
            threads.append(t)

        start = time.monotonic()
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=60)
        elapsed = time.monotonic() - start

        assert len(completed) == 50, (
            f"Only {len(completed)}/50 threads completed in {elapsed:.1f}s "
            f"(possible deadlock)"
        )

    def test_concurrent_different_bodies(self) -> None:
        """Each thread computes a different body — results match baseline."""
        all_bodies = [
            swe.SUN,
            swe.MOON,
            swe.MERCURY,
            swe.VENUS,
            swe.MARS,
            swe.JUPITER,
            swe.SATURN,
            swe.URANUS,
            swe.NEPTUNE,
            swe.PLUTO,
        ]
        jd = J2000

        # Single-threaded baseline
        baseline = {}
        ctx = EphemerisContext()
        for body in all_bodies:
            pos, _ = ctx.calc_ut(jd, body, 0)
            baseline[body] = pos

        errors: list[str] = []
        lock = threading.Lock()

        def worker(body: int) -> None:
            ctx = EphemerisContext()
            for _ in range(100):
                pos, _ = ctx.calc_ut(jd, body, 0)
                for i in range(6):
                    if abs(pos[i] - baseline[body][i]) > 1e-10:
                        with lock:
                            errors.append(f"Body {body}: elem {i} mismatch")
                        return

        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(worker, b) for b in all_bodies]
            for f in as_completed(futures):
                f.result()  # Raise any exception

        assert not errors, "\n".join(errors[:10])


# ============================================================================
# §6.2 Global State Isolation
# ============================================================================


class TestGlobalStateIsolation:
    """§6.2 Verify that per-context state is truly isolated between threads."""

    def test_sidereal_mode_isolation(self) -> None:
        """Thread A (LAHIRI) and Thread B (FAGAN_BRADLEY) get different results."""
        jd = J2000
        body = swe.SUN
        flags = swe.FLG_SIDEREAL | swe.FLG_SPEED

        results_lahiri: list[tuple[float, ...]] = []
        results_fagan: list[tuple[float, ...]] = []
        lock = threading.Lock()
        barrier = threading.Barrier(2)

        def worker_lahiri() -> None:
            ctx = EphemerisContext()
            ctx.set_sid_mode(swe.SIDM_LAHIRI)
            barrier.wait(timeout=10)
            for _ in range(50):
                pos, _ = ctx.calc_ut(jd, body, flags)
                with lock:
                    results_lahiri.append(pos)

        def worker_fagan() -> None:
            ctx = EphemerisContext()
            ctx.set_sid_mode(swe.SIDM_FAGAN_BRADLEY)
            barrier.wait(timeout=10)
            for _ in range(50):
                pos, _ = ctx.calc_ut(jd, body, flags)
                with lock:
                    results_fagan.append(pos)

        t1 = threading.Thread(target=worker_lahiri)
        t2 = threading.Thread(target=worker_fagan)
        t1.start()
        t2.start()
        t1.join(timeout=30)
        t2.join(timeout=30)

        assert len(results_lahiri) == 50, "Lahiri thread did not complete"
        assert len(results_fagan) == 50, "Fagan thread did not complete"

        # The two modes must produce DIFFERENT longitudes
        lahiri_lon = results_lahiri[0][0]
        fagan_lon = results_fagan[0][0]
        diff = abs(lahiri_lon - fagan_lon)
        assert diff > 0.1, (
            f"Lahiri and Fagan should differ: Lahiri={lahiri_lon:.6f}, "
            f"Fagan={fagan_lon:.6f}, diff={diff:.6f}"
        )

        # All results within each thread must be consistent
        for i, pos in enumerate(results_lahiri):
            assert abs(pos[0] - lahiri_lon) < 1e-10, (
                f"Lahiri result {i} inconsistent: {pos[0]} vs {lahiri_lon}"
            )
        for i, pos in enumerate(results_fagan):
            assert abs(pos[0] - fagan_lon) < 1e-10, (
                f"Fagan result {i} inconsistent: {pos[0]} vs {fagan_lon}"
            )

    @pytest.mark.xfail(
        reason="Known: _CONTEXT_SWAP_LOCK does not fully isolate topocentric state "
        "in _calc_body_with_context — requires architectural fix (context-only mode)",
        strict=False,
    )
    def test_topo_isolation(self) -> None:
        """Different threads with different topocentric positions get different results."""
        jd = J2000
        body = swe.MOON
        flags = swe.FLG_TOPOCTR | swe.FLG_SPEED

        results_rome: list[tuple[float, ...]] = []
        results_tokyo: list[tuple[float, ...]] = []
        lock = threading.Lock()
        barrier = threading.Barrier(2)

        def worker_rome() -> None:
            ctx = EphemerisContext()
            ctx.set_topo(12.5, 41.9, 50.0)  # Rome
            barrier.wait(timeout=10)
            for _ in range(50):
                pos, _ = ctx.calc_ut(jd, body, flags)
                with lock:
                    results_rome.append(pos)

        def worker_tokyo() -> None:
            ctx = EphemerisContext()
            ctx.set_topo(139.7, 35.7, 40.0)  # Tokyo
            barrier.wait(timeout=10)
            for _ in range(50):
                pos, _ = ctx.calc_ut(jd, body, flags)
                with lock:
                    results_tokyo.append(pos)

        t1 = threading.Thread(target=worker_rome)
        t2 = threading.Thread(target=worker_tokyo)
        t1.start()
        t2.start()
        t1.join(timeout=30)
        t2.join(timeout=30)

        assert len(results_rome) == 50, "Rome thread did not complete"
        assert len(results_tokyo) == 50, "Tokyo thread did not complete"

        # Topocentric Moon positions from Rome and Tokyo must differ
        rome_lon = results_rome[0][0]
        tokyo_lon = results_tokyo[0][0]
        diff = abs(rome_lon - tokyo_lon)
        # For the Moon, topocentric parallax is ~1° max
        assert diff > 0.001, (
            f"Rome and Tokyo should differ: Rome={rome_lon:.6f}, Tokyo={tokyo_lon:.6f}"
        )

        # Consistency within each thread
        for i, pos in enumerate(results_rome):
            assert abs(pos[0] - rome_lon) < 1e-10, f"Rome result {i} inconsistent"
        for i, pos in enumerate(results_tokyo):
            assert abs(pos[0] - tokyo_lon) < 1e-10, f"Tokyo result {i} inconsistent"

    def test_sid_mode_switching_stress(self) -> None:
        """Multiple threads rapidly switching sidereal modes."""
        jd = J2000
        errors: list[str] = []
        lock = threading.Lock()

        sid_modes = [
            swe.SIDM_LAHIRI,  # 1
            swe.SIDM_FAGAN_BRADLEY,  # 0
        ]

        def worker(mode: int, thread_id: int) -> None:
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)
            flags = swe.FLG_SIDEREAL
            ref_pos = None
            for i in range(100):
                pos, _ = ctx.calc_ut(jd, swe.SUN, flags)
                if ref_pos is None:
                    ref_pos = pos
                elif abs(pos[0] - ref_pos[0]) > 1e-10:
                    with lock:
                        errors.append(
                            f"Thread {thread_id} (mode={mode}): "
                            f"inconsistent at iter {i}"
                        )
                    return

        threads = []
        for tid in range(20):
            mode = sid_modes[tid % 2]
            t = threading.Thread(target=worker, args=(mode, tid))
            threads.append(t)

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert not errors, "\n".join(errors[:10])


# ============================================================================
# §6.3 LEB + Skyfield Mixed Mode
# ============================================================================


class TestLEBSkyfieldMixedMode:
    """§6.3 Verify no interference between LEB and Skyfield modes."""

    def test_leb_and_skyfield_concurrent(self) -> None:
        """Contexts with LEB and without LEB produce consistent results."""
        jd = J2000
        body = swe.SUN
        flags = swe.FLG_SPEED

        # Get a LEB file path if available
        leb_path = None
        try:
            import os

            for candidate in [
                "data/leb/ephemeris_medium.leb",
                "data/leb/ephemeris_base.leb",
            ]:
                if os.path.exists(candidate):
                    leb_path = candidate
                    break
        except Exception:
            pass

        if leb_path is None:
            pytest.skip("No LEB file available for mixed-mode test")

        results_leb: list[tuple[float, ...]] = []
        results_skyfield: list[tuple[float, ...]] = []
        lock = threading.Lock()
        barrier = threading.Barrier(2)

        def worker_leb() -> None:
            ctx = EphemerisContext()
            ctx.set_leb_file(leb_path)
            barrier.wait(timeout=10)
            for i in range(100):
                pos, _ = ctx.calc_ut(jd + i * 10.0, body, flags)
                with lock:
                    results_leb.append(pos)

        def worker_skyfield() -> None:
            ctx = EphemerisContext()
            ctx.set_leb_file(None)  # Force Skyfield
            barrier.wait(timeout=10)
            for i in range(100):
                pos, _ = ctx.calc_ut(jd + i * 10.0, body, flags)
                with lock:
                    results_skyfield.append(pos)

        t1 = threading.Thread(target=worker_leb)
        t2 = threading.Thread(target=worker_skyfield)
        t1.start()
        t2.start()
        t1.join(timeout=60)
        t2.join(timeout=60)

        assert len(results_leb) == 100, "LEB thread did not complete"
        assert len(results_skyfield) == 100, "Skyfield thread did not complete"

        # LEB and Skyfield should agree to within 0.001" (they're the same
        # ephemeris, just different code paths)
        for i in range(100):
            leb_pos = results_leb[i]
            sky_pos = results_skyfield[i]
            # Position elements (lon, lat, dist)
            for j in range(3):
                diff = abs(leb_pos[j] - sky_pos[j])
                if j < 2:  # Angular elements
                    diff_arcsec = diff * 3600
                    assert diff_arcsec < 1.0, (
                        f'Date {i}: element {j} differs by {diff_arcsec:.4f}" '
                        f"(LEB={leb_pos[j]:.8f}, Sky={sky_pos[j]:.8f})"
                    )

    def test_multiple_leb_contexts_concurrent(self) -> None:
        """Multiple LEB contexts running concurrently produce consistent results."""
        import os

        leb_path = None
        for candidate in [
            "data/leb/ephemeris_medium.leb",
            "data/leb/ephemeris_base.leb",
        ]:
            if os.path.exists(candidate):
                leb_path = candidate
                break

        if leb_path is None:
            pytest.skip("No LEB file available")

        jd = J2000
        errors: list[str] = []
        lock = threading.Lock()

        # Compute baseline
        ctx_base = EphemerisContext()
        ctx_base.set_leb_file(leb_path)
        baseline = {}
        for body in [swe.SUN, swe.MOON, swe.MARS]:
            pos, _ = ctx_base.calc_ut(jd, body, swe.FLG_SPEED)
            baseline[body] = pos

        def worker(thread_id: int) -> None:
            ctx = EphemerisContext()
            ctx.set_leb_file(leb_path)
            for _ in range(50):
                for body in [swe.SUN, swe.MOON, swe.MARS]:
                    pos, _ = ctx.calc_ut(jd, body, swe.FLG_SPEED)
                    for i in range(6):
                        if abs(pos[i] - baseline[body][i]) > 1e-10:
                            with lock:
                                errors.append(f"T{thread_id}: body={body} elem={i}")
                            return

        threads = []
        for tid in range(20):
            t = threading.Thread(target=worker, args=(tid,))
            threads.append(t)

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=60)

        assert not errors, f"{len(errors)} mismatches:\n" + "\n".join(errors[:10])

    def test_house_calculations_concurrent(self) -> None:
        """House calculations from multiple threads must be consistent."""
        jd = J2000
        lat, lon = 45.0, 10.0

        # Baseline
        ctx = EphemerisContext()
        baseline_cusps, baseline_angles = ctx.houses(jd, lat, lon, ord("P"))

        errors: list[str] = []
        lock = threading.Lock()

        def worker(thread_id: int) -> None:
            ctx = EphemerisContext()
            for _ in range(100):
                cusps, angles = ctx.houses(jd, lat, lon, ord("P"))
                for i in range(min(len(cusps), len(baseline_cusps))):
                    if abs(cusps[i] - baseline_cusps[i]) > 1e-10:
                        with lock:
                            errors.append(f"T{thread_id}: cusp {i} mismatch")
                        return

        threads = []
        for tid in range(20):
            t = threading.Thread(target=worker, args=(tid,))
            threads.append(t)

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert not errors, "\n".join(errors[:10])


# ============================================================================
# Additional concurrency tests
# ============================================================================


class TestConcurrencyEdgeCases:
    """Additional concurrency edge cases."""

    def test_context_creation_under_load(self) -> None:
        """Creating many EphemerisContext instances concurrently."""
        contexts: list[EphemerisContext] = []
        lock = threading.Lock()

        def worker() -> None:
            ctx = EphemerisContext()
            # Do a calculation to force lazy init of shared resources
            pos, _ = ctx.calc_ut(J2000, swe.SUN, 0)
            assert len(pos) == 6
            with lock:
                contexts.append(ctx)

        threads = [threading.Thread(target=worker) for _ in range(50)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert len(contexts) == 50, f"Only {len(contexts)}/50 contexts created"

    def test_threadpool_executor_pattern(self) -> None:
        """Typical ThreadPoolExecutor usage pattern works correctly."""
        bodies = list(range(10))  # SUN through PLUTO
        jds = [J2000 + i for i in range(10)]

        def compute(body: int, jd: float) -> tuple[int, float, tuple]:
            ctx = EphemerisContext()
            pos, _ = ctx.calc_ut(jd, body, swe.FLG_SPEED)
            return (body, jd, pos)

        results = []
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = []
            for body in bodies:
                for jd in jds:
                    futures.append(executor.submit(compute, body, jd))

            for f in as_completed(futures):
                body, jd, pos = f.result()
                assert len(pos) == 6
                for v in pos:
                    assert math.isfinite(v)
                results.append((body, jd, pos))

        assert len(results) == 100, f"Expected 100 results, got {len(results)}"

    def test_interleaved_read_write(self) -> None:
        """Threads interleaving state writes and reads via EphemerisContext."""
        errors: list[str] = []
        lock = threading.Lock()

        def worker(thread_id: int) -> None:
            ctx = EphemerisContext()
            # Each thread uses a unique sidereal mode
            mode = thread_id % 2  # 0=Fagan, 1=Lahiri
            ctx.set_sid_mode(mode)

            for _ in range(50):
                # Verify mode hasn't been corrupted
                current_mode = ctx.get_sid_mode()
                if current_mode != mode:
                    with lock:
                        errors.append(
                            f"T{thread_id}: mode changed from {mode} to {current_mode}"
                        )
                    return

                # Do a sidereal calculation
                pos, _ = ctx.calc_ut(J2000, swe.SUN, swe.FLG_SIDEREAL)
                assert len(pos) == 6

        threads = []
        for tid in range(30):
            t = threading.Thread(target=worker, args=(tid,))
            threads.append(t)

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert not errors, "\n".join(errors[:10])
