"""
Advanced thread safety verification tests for EphemerisContext.

This module provides rigorous tests that verify EphemerisContext provides true
thread safety for multi-threaded applications. These tests go beyond basic
concurrent access testing to verify:

1. Lock correctness - The context swap lock properly serializes access
2. No deadlocks - High contention doesn't cause deadlocks
3. State isolation - Each context's state is truly isolated
4. Result consistency - Concurrent calculations produce correct results
5. Memory visibility - State changes are visible across threads
6. Stress resilience - Thread safety holds under extreme load

These tests are designed to catch race conditions that may only manifest
under specific timing conditions.
"""

import concurrent.futures
import random
import threading
import time
from typing import Dict, List, Tuple

import pytest

from libephemeris import EphemerisContext
from libephemeris.constants import (
    JUPITER,
    MERCURY,
    MOON,
    SATURN,
    SUN,
    VENUS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_TOPOCTR,
    SIDM_FAGAN_BRADLEY,
    SIDM_KRISHNAMURTI,
    SIDM_LAHIRI,
    SIDM_RAMAN,
    SIDM_TRUE_CITRA,
    SIDM_YUKTESHWAR,
)

pytestmark = pytest.mark.slow


# =============================================================================
# TEST: Lock Correctness Verification
# =============================================================================


class TestLockCorrectness:
    """
    Verify that _CONTEXT_SWAP_LOCK correctly serializes access.

    These tests ensure the RLock properly protects the save-set-restore cycle
    and prevents interleaving of operations.
    """

    def test_lock_exists_and_is_reentrant(self):
        """Verify the context swap lock is an RLock (reentrant lock)."""
        from libephemeris import state

        assert hasattr(state, "_CONTEXT_SWAP_LOCK")

        # Verify it's reentrant by acquiring it twice from the same thread
        lock = state._CONTEXT_SWAP_LOCK
        acquired_first = lock.acquire(blocking=False)
        assert acquired_first, "Should acquire lock first time"

        # RLock should allow re-acquisition from same thread
        acquired_second = lock.acquire(blocking=False)
        assert acquired_second, "RLock should allow reentrant acquisition"

        # Release both acquisitions
        lock.release()
        lock.release()

    def test_shared_lock_exists_and_is_reentrant(self):
        """Verify the shared resource lock in context.py is an RLock."""
        from libephemeris.context import _SHARED_LOCK

        # Verify it's reentrant
        acquired_first = _SHARED_LOCK.acquire(blocking=False)
        assert acquired_first, "Should acquire lock first time"

        acquired_second = _SHARED_LOCK.acquire(blocking=False)
        assert acquired_second, "RLock should allow reentrant acquisition"

        _SHARED_LOCK.release()
        _SHARED_LOCK.release()

    def test_context_swap_is_atomic(self):
        """
        Verify that context swap operations are atomic.

        This test creates a scenario where multiple threads attempt to
        interleave their context swaps. If the lock is working correctly,
        each thread should see consistent state throughout its calculation.
        """

        jd = 2451545.0
        num_iterations = 100
        errors: List[str] = []
        errors_lock = threading.Lock()

        # Two very different ayanamsha configurations
        configs = [
            (SIDM_LAHIRI, "Lahiri"),
            (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
            (SIDM_RAMAN, "Raman"),
            (SIDM_KRISHNAMURTI, "Krishnamurti"),
        ]

        # Pre-calculate expected values for each config
        expected: Dict[int, float] = {}
        for mode, _ in configs:
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)
            pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
            expected[mode] = pos[0]

        def worker(thread_id: int, mode: int, mode_name: str):
            """Each worker calculates with its own mode and verifies correctness."""
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(mode)

                for i in range(num_iterations):
                    pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)

                    # Verify the result matches what we expect for this mode
                    if abs(pos[0] - expected[mode]) > 1e-6:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} ({mode_name}) iteration {i}: "
                                f"got {pos[0]:.6f}, expected {expected[mode]:.6f}"
                            )
                        return
            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id} ({mode_name}): {e}")

        # Launch threads with different configurations simultaneously
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(configs)) as ex:
            futures = [
                ex.submit(worker, i, mode, name)
                for i, (mode, name) in enumerate(configs)
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Atomic swap violations detected: {errors}"

    def test_nested_context_calculations(self):
        """
        Test that nested context operations don't deadlock.

        Since we use RLock, the same thread should be able to acquire the lock
        multiple times during nested calculations.
        """
        jd = 2451545.0
        errors: List[str] = []

        # This simulates a complex calculation that might internally call
        # other context-aware functions
        def nested_calculation():
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(SIDM_LAHIRI)

                # First calculation
                pos1, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)

                # Second calculation (same thread, same lock)
                pos2, _ = ctx.calc_ut(jd, MOON, FLG_SIDEREAL)

                # House calculation (also uses the lock)
                cusps, ascmc = ctx.houses(jd, 41.9, 12.5, ord("P"))

                # Planet-centric calculation (also uses the lock)
                # Note: Using SUN as center since MARS may not be in all kernels
                pos3, _ = ctx.calc_pctr(jd, MOON, SUN, FLG_SIDEREAL)

                return True
            except Exception as e:
                errors.append(f"Nested calculation error: {e}")
                return False

        # Run multiple nested calculations concurrently
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
            futures = [ex.submit(nested_calculation) for _ in range(16)]
            results = [f.result() for f in futures]

        assert not errors, f"Nested calculation errors: {errors}"
        assert all(results), "All nested calculations should succeed"


# =============================================================================
# TEST: No Deadlocks Under Contention
# =============================================================================


class TestNoDeadlocks:
    """
    Verify that high contention doesn't cause deadlocks.

    These tests create scenarios with many threads competing for the lock
    to verify the implementation doesn't deadlock.
    """

    def test_high_contention_no_deadlock(self):
        """
        Test that 100 concurrent threads don't deadlock.

        Each thread performs calculations while other threads compete for
        the same lock. All threads should eventually complete.
        """
        num_threads = 100
        jd = 2451545.0
        timeout_seconds = 60  # If it takes longer, likely deadlocked
        completed = [0]
        completed_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                mode = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN][
                    thread_id % 3
                ]
                ctx.set_sid_mode(mode)

                # Perform a few calculations
                for _ in range(10):
                    ctx.calc_ut(jd + thread_id * 0.1, SUN, FLG_SIDEREAL)
                    ctx.calc_ut(jd + thread_id * 0.1, MOON, FLG_SIDEREAL)

                with completed_lock:
                    completed[0] += 1

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        start_time = time.time()

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            # Wait with timeout to detect deadlocks
            done, not_done = concurrent.futures.wait(futures, timeout=timeout_seconds)

        elapsed = time.time() - start_time

        if not_done:
            pytest.fail(
                f"Potential deadlock: {len(not_done)} threads did not complete "
                f"within {timeout_seconds}s"
            )

        assert not errors, f"Errors occurred: {errors[:10]}"
        assert completed[0] == num_threads, (
            f"Expected {num_threads} completions, got {completed[0]}"
        )
        assert elapsed < timeout_seconds, f"Took too long: {elapsed:.1f}s"

    def test_mixed_operation_contention(self):
        """
        Test mixed read/write operations don't deadlock.

        Different threads perform different types of operations that all
        require the context swap lock.
        """
        num_threads = 50
        jd = 2451545.0
        timeout_seconds = 30
        completed = [0]
        completed_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        def planet_worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(SIDM_LAHIRI)
                for _ in range(5):
                    ctx.calc_ut(jd, SUN, FLG_SIDEREAL | FLG_SPEED)
                with completed_lock:
                    completed[0] += 1
            except Exception as e:
                with errors_lock:
                    errors.append(f"Planet worker {thread_id}: {e}")

        def house_worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(SIDM_FAGAN_BRADLEY)
                for _ in range(5):
                    ctx.houses(jd, 41.9, 12.5, ord("P"))
                with completed_lock:
                    completed[0] += 1
            except Exception as e:
                with errors_lock:
                    errors.append(f"House worker {thread_id}: {e}")

        def pctr_worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(SIDM_RAMAN)
                for _ in range(5):
                    # Use SUN as center since MARS may not be in all kernels
                    ctx.calc_pctr(jd, MOON, SUN, FLG_SIDEREAL)
                with completed_lock:
                    completed[0] += 1
            except Exception as e:
                with errors_lock:
                    errors.append(f"Pctr worker {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = []
            for i in range(num_threads):
                if i % 3 == 0:
                    futures.append(ex.submit(planet_worker, i))
                elif i % 3 == 1:
                    futures.append(ex.submit(house_worker, i))
                else:
                    futures.append(ex.submit(pctr_worker, i))

            done, not_done = concurrent.futures.wait(futures, timeout=timeout_seconds)

        if not_done:
            pytest.fail(f"Potential deadlock: {len(not_done)} threads did not complete")

        assert not errors, f"Errors occurred: {errors}"
        assert completed[0] == num_threads


# =============================================================================
# TEST: State Isolation Verification
# =============================================================================


class TestStateIsolationVerification:
    """
    Rigorously verify that each EphemerisContext has truly isolated state.

    These tests verify that modifications to one context never affect another.
    """

    def test_sidereal_mode_isolation_under_stress(self):
        """
        Verify sidereal mode isolation with many concurrent modifications.

        Each thread creates a context, sets a unique mode, performs
        calculations, and verifies its mode wasn't changed by other threads.
        """
        num_threads = 20
        num_iterations = 50
        errors: List[str] = []
        errors_lock = threading.Lock()

        modes = [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
            SIDM_KRISHNAMURTI,
            SIDM_YUKTESHWAR,
            SIDM_TRUE_CITRA,
        ]

        def worker(thread_id: int):
            mode = modes[thread_id % len(modes)]
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)

            try:
                for i in range(num_iterations):
                    # Verify mode before calculation
                    if ctx.get_sid_mode() != mode:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} iter {i}: mode changed "
                                f"before calc (expected {mode}, got {ctx.get_sid_mode()})"
                            )
                        return

                    # Perform calculation
                    ctx.calc_ut(2451545.0, SUN, FLG_SIDEREAL)

                    # Verify mode after calculation
                    if ctx.get_sid_mode() != mode:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} iter {i}: mode changed "
                                f"after calc (expected {mode}, got {ctx.get_sid_mode()})"
                            )
                        return

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"State isolation violations: {errors}"

    def test_topo_isolation_under_stress(self):
        """
        Verify topocentric location isolation with concurrent modifications.

        Each thread sets a unique location and verifies it isn't affected
        by other threads setting different locations.
        """
        locations = [
            (0.0, 0.0, 0, "Equator"),
            (12.5, 41.9, 50, "Rome"),
            (-0.1, 51.5, 10, "London"),
            (139.7, 35.7, 40, "Tokyo"),
            (-74.0, 40.7, 10, "New York"),
            (151.2, -33.9, 5, "Sydney"),
        ]
        num_threads = len(locations)
        num_iterations = 50
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            lon, lat, alt, name = locations[thread_id % len(locations)]
            ctx = EphemerisContext()
            ctx.set_topo(lon, lat, alt)

            try:
                for i in range(num_iterations):
                    # Verify topo before calculation
                    topo = ctx.get_topo()
                    if topo is None:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} ({name}) iter {i}: topo is None"
                            )
                        return

                    if abs(topo.latitude.degrees - lat) > 0.001:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} ({name}) iter {i}: lat changed "
                                f"(expected {lat}, got {topo.latitude.degrees})"
                            )
                        return

                    # Perform calculation
                    ctx.calc_ut(2451545.0, SUN, FLG_TOPOCTR | FLG_SPEED)

                    # Verify topo after calculation
                    topo = ctx.get_topo()
                    if abs(topo.latitude.degrees - lat) > 0.001:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} ({name}) iter {i}: lat changed "
                                f"after calc (expected {lat}, got {topo.latitude.degrees})"
                            )
                        return

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id} ({name}): {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Topo isolation violations: {errors}"

    def test_angles_cache_isolation(self):
        """
        Verify that angles cache is isolated per context.

        Each thread sets different cached angles and verifies they don't
        leak to other contexts.
        """
        num_threads = 10
        num_iterations = 50
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            ctx = EphemerisContext()
            # Set unique angles for this thread
            angles = {
                "Sun": float(thread_id * 10),
                "Moon": float(thread_id * 20),
                "Asc": float(thread_id * 30),
            }
            ctx.set_angles_cache(angles)

            try:
                for i in range(num_iterations):
                    # Verify our angles haven't been modified
                    cached = ctx.get_angles_cache()
                    if cached.get("Sun") != angles["Sun"]:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} iter {i}: Sun angle changed "
                                f"(expected {angles['Sun']}, got {cached.get('Sun')})"
                            )
                        return

                    # Do a calculation that uses global state
                    ctx.calc_ut(2451545.0, SUN, FLG_SIDEREAL)

                    # Verify our cache is still intact
                    cached = ctx.get_angles_cache()
                    if cached.get("Moon") != angles["Moon"]:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} iter {i}: Moon angle changed "
                                f"after calc"
                            )
                        return

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Angles cache isolation violations: {errors}"


# =============================================================================
# TEST: Result Consistency Under Concurrency
# =============================================================================


class TestResultConsistencyUnderConcurrency:
    """
    Verify that concurrent calculations produce mathematically correct results.

    These tests compare concurrent results against known-good sequential
    calculations to detect any corruption due to race conditions.
    """

    def test_sidereal_positions_match_sequential(self):
        """
        Verify concurrent sidereal calculations match sequential ones.

        Pre-calculate expected results sequentially, then verify concurrent
        calculations produce identical results.
        """
        jd = 2451545.0
        planets = [
            SUN,
            MOON,
            MERCURY,
            VENUS,
            # Note: MARS may not be directly available in all ephemeris kernels
            JUPITER,
            SATURN,
        ]
        modes = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN]

        # Pre-calculate expected results
        expected: Dict[Tuple[int, int], float] = {}
        for mode in modes:
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)
            for planet in planets:
                pos, _ = ctx.calc_ut(jd, planet, FLG_SIDEREAL)
                expected[(mode, planet)] = pos[0]

        num_iterations = 100
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            mode = modes[thread_id % len(modes)]
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)

            try:
                for i in range(num_iterations):
                    planet = planets[i % len(planets)]
                    pos, _ = ctx.calc_ut(jd, planet, FLG_SIDEREAL)

                    exp = expected[(mode, planet)]
                    if abs(pos[0] - exp) > 1e-8:
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} iter {i}: mode={mode} planet={planet} "
                                f"got {pos[0]:.10f}, expected {exp:.10f}"
                            )
                        return

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            futures = [ex.submit(worker, i) for i in range(12)]
            concurrent.futures.wait(futures)

        assert not errors, f"Result consistency errors: {errors}"

    def test_house_cusps_match_sequential(self):
        """
        Verify concurrent house calculations match sequential ones.
        """
        jd = 2451545.0
        locations = [
            (12.5, 41.9),  # Rome
            (-0.1, 51.5),  # London
            (139.7, 35.7),  # Tokyo
        ]
        modes = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY]

        # Pre-calculate expected results
        expected: Dict[Tuple[int, float, float], Tuple[float, ...]] = {}
        for mode in modes:
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)
            for lon, lat in locations:
                cusps, ascmc = ctx.houses(jd, lat, lon, ord("P"))
                expected[(mode, lon, lat)] = cusps

        num_iterations = 50
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            mode = modes[thread_id % len(modes)]
            lon, lat = locations[thread_id % len(locations)]
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)

            try:
                for i in range(num_iterations):
                    cusps, _ = ctx.houses(jd, lat, lon, ord("P"))
                    exp_cusps = expected[(mode, lon, lat)]

                    for j, (got, exp) in enumerate(zip(cusps, exp_cusps)):
                        if abs(got - exp) > 1e-6:
                            with errors_lock:
                                errors.append(
                                    f"Thread {thread_id} iter {i}: cusp {j} "
                                    f"got {got:.6f}, expected {exp:.6f}"
                                )
                            return

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
            futures = [ex.submit(worker, i) for i in range(6)]
            concurrent.futures.wait(futures)

        assert not errors, f"House calculation errors: {errors}"

    def test_repeated_same_calculation_consistent(self):
        """
        Verify the same calculation repeated by multiple threads gives
        identical results every time.
        """
        jd = 2451545.0
        num_threads = 20
        num_iterations = 100
        results: List[float] = []
        results_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        # All threads use the same configuration
        mode = SIDM_LAHIRI
        planet = SUN

        def worker(thread_id: int):
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)

            try:
                for _ in range(num_iterations):
                    pos, _ = ctx.calc_ut(jd, planet, FLG_SIDEREAL)
                    with results_lock:
                        results.append(pos[0])

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"

        # All results should be identical
        if results:
            reference = results[0]
            for i, r in enumerate(results):
                if abs(r - reference) > 1e-10:
                    pytest.fail(
                        f"Result {i} differs: {r:.10f} vs reference {reference:.10f}"
                    )


# =============================================================================
# TEST: Memory Visibility
# =============================================================================


class TestMemoryVisibility:
    """
    Verify that state changes are properly visible across threads.

    Python's GIL provides some memory ordering guarantees, but these tests
    verify our implementation doesn't have any subtle visibility issues.
    """

    def test_context_configuration_visible_in_calculation(self):
        """
        Verify that configuration set before calculation is correctly used.

        This tests the fundamental contract that settings applied to a context
        are visible during the calculation.
        """
        jd = 2451545.0
        num_threads = 10
        errors: List[str] = []
        errors_lock = threading.Lock()

        # Pre-calculate expected difference between two modes
        ctx_lahiri = EphemerisContext()
        ctx_lahiri.set_sid_mode(SIDM_LAHIRI)
        pos_lahiri, _ = ctx_lahiri.calc_ut(jd, SUN, FLG_SIDEREAL)

        ctx_fagan = EphemerisContext()
        ctx_fagan.set_sid_mode(SIDM_FAGAN_BRADLEY)
        pos_fagan, _ = ctx_fagan.calc_ut(jd, SUN, FLG_SIDEREAL)

        # The two should be different
        assert abs(pos_lahiri[0] - pos_fagan[0]) > 0.1, (
            "Modes should produce different results"
        )

        def worker(thread_id: int):
            # Alternate between modes
            use_lahiri = thread_id % 2 == 0
            expected = pos_lahiri[0] if use_lahiri else pos_fagan[0]

            try:
                ctx = EphemerisContext()
                if use_lahiri:
                    ctx.set_sid_mode(SIDM_LAHIRI)
                else:
                    ctx.set_sid_mode(SIDM_FAGAN_BRADLEY)

                # Verify the calculation uses the mode we just set
                for i in range(100):
                    pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
                    if abs(pos[0] - expected) > 1e-6:
                        with errors_lock:
                            mode_name = "Lahiri" if use_lahiri else "Fagan-Bradley"
                            errors.append(
                                f"Thread {thread_id} ({mode_name}) iter {i}: "
                                f"got {pos[0]:.6f}, expected {expected:.6f}"
                            )
                        return

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Memory visibility errors: {errors}"


# =============================================================================
# TEST: Stress Tests
# =============================================================================


class TestStressConditions:
    """
    Extreme stress tests to verify thread safety holds under pressure.
    """

    def test_burst_concurrent_calculations(self):
        """
        Test a burst of many calculations happening simultaneously.

        This creates maximum contention on the context swap lock.
        """
        num_threads = 100
        calculations_per_thread = 10
        jd = 2451545.0
        completed_calculations = [0]
        count_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        modes = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN]
        # Note: Avoid MARS which may not be directly available in all kernels
        planets = [SUN, MOON, VENUS, JUPITER]

        def worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(modes[thread_id % len(modes)])

                for i in range(calculations_per_thread):
                    planet = planets[i % len(planets)]
                    pos, _ = ctx.calc_ut(jd + i * 0.1, planet, FLG_SIDEREAL)

                    if not (0 <= pos[0] < 360):
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} calc {i}: invalid position {pos[0]}"
                            )
                        return

                with count_lock:
                    completed_calculations[0] += calculations_per_thread

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        start = time.time()

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        elapsed = time.time() - start
        expected = num_threads * calculations_per_thread

        assert not errors, f"Errors: {errors[:20]}"
        assert completed_calculations[0] == expected, (
            f"Expected {expected} calculations, got {completed_calculations[0]}"
        )
        # Should complete in reasonable time (not deadlocked)
        assert elapsed < 120, f"Burst test took too long: {elapsed:.1f}s"

    def test_randomized_operations(self):
        """
        Test with randomized operations to catch edge cases.

        Each thread randomly chooses operations and configurations to
        maximize the chance of finding race conditions.
        """
        num_threads = 30
        operations_per_thread = 50
        jd = 2451545.0
        errors: List[str] = []
        errors_lock = threading.Lock()
        completed = [0]
        completed_lock = threading.Lock()

        modes = [
            SIDM_LAHIRI,
            SIDM_FAGAN_BRADLEY,
            SIDM_RAMAN,
            SIDM_KRISHNAMURTI,
        ]
        # Note: Avoid MARS which may not be directly available in all kernels
        planets = [SUN, MOON, MERCURY, VENUS, JUPITER]
        locations = [
            (12.5, 41.9, 50),
            (-0.1, 51.5, 10),
            (139.7, 35.7, 40),
        ]

        def worker(thread_id: int):
            # Each thread gets its own random seed for reproducibility
            rng = random.Random(thread_id)
            ctx = EphemerisContext()

            try:
                for _ in range(operations_per_thread):
                    operation = rng.choice(["calc", "calc_sidereal", "houses", "pctr"])

                    if operation == "calc":
                        planet = rng.choice(planets)
                        ctx.calc_ut(jd + rng.random() * 100, planet, FLG_SPEED)

                    elif operation == "calc_sidereal":
                        ctx.set_sid_mode(rng.choice(modes))
                        planet = rng.choice(planets)
                        ctx.calc_ut(jd + rng.random() * 100, planet, FLG_SIDEREAL)

                    elif operation == "houses":
                        ctx.set_sid_mode(rng.choice(modes))
                        lon, lat, _ = rng.choice(locations)
                        ctx.houses(jd + rng.random() * 100, lat, lon, ord("P"))

                    elif operation == "pctr":
                        ctx.set_sid_mode(rng.choice(modes))
                        ctx.calc_pctr(
                            jd + rng.random() * 100,
                            rng.choice([MOON, VENUS]),
                            rng.choice([SUN, JUPITER]),
                            FLG_SIDEREAL,
                        )

                with completed_lock:
                    completed[0] += 1

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors in randomized test: {errors}"
        assert completed[0] == num_threads, (
            f"Expected {num_threads} threads to complete, got {completed[0]}"
        )

    def test_rapid_context_creation_and_use(self):
        """
        Test rapid creation and use of many contexts.

        This stresses both the shared resource initialization and the
        context swap mechanism.
        """
        num_threads = 20
        contexts_per_thread = 100
        jd = 2451545.0
        errors: List[str] = []
        errors_lock = threading.Lock()
        context_count = [0]
        count_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                mode = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY][thread_id % 2]

                for i in range(contexts_per_thread):
                    # Create a fresh context each iteration
                    ctx = EphemerisContext()
                    ctx.set_sid_mode(mode)

                    # Use it for a calculation
                    pos, _ = ctx.calc_ut(jd + i * 0.01, SUN, FLG_SIDEREAL)

                    if not (0 <= pos[0] < 360):
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id} ctx {i}: invalid position"
                            )
                        return

                with count_lock:
                    context_count[0] += contexts_per_thread

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        expected = num_threads * contexts_per_thread
        assert not errors, f"Errors: {errors}"
        assert context_count[0] == expected, (
            f"Expected {expected} contexts, got {context_count[0]}"
        )


# =============================================================================
# TEST: Shared Resource Thread Safety
# =============================================================================


class TestSharedResourceThreadSafety:
    """
    Verify thread-safe initialization and access to shared resources.

    These tests verify the double-checked locking pattern used for
    lazy initialization of shared resources.
    """

    def test_shared_loader_initialization(self):
        """
        Verify multiple threads can safely initialize the shared loader.
        """

        # Reset shared resources
        EphemerisContext.close()

        num_threads = 20
        loaders: List[object] = []
        loaders_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                loader = ctx.get_loader()
                with loaders_lock:
                    loaders.append(id(loader))
            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        assert len(loaders) == num_threads

        # All should get the same shared loader instance
        unique_loaders = set(loaders)
        assert len(unique_loaders) == 1, (
            f"Expected 1 shared loader, got {len(unique_loaders)} different instances"
        )

    def test_shared_timescale_initialization(self):
        """
        Verify multiple threads can safely initialize the shared timescale.
        """
        # Reset shared resources
        EphemerisContext.close()

        num_threads = 20
        timescales: List[int] = []
        timescales_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                ts = ctx.get_timescale()
                with timescales_lock:
                    timescales.append(id(ts))
            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        assert len(timescales) == num_threads

        # All should get the same shared timescale instance
        unique_timescales = set(timescales)
        assert len(unique_timescales) == 1, (
            f"Expected 1 shared timescale, got {len(unique_timescales)} instances"
        )

    def test_shared_planets_initialization(self):
        """
        Verify multiple threads can safely initialize the shared ephemeris.
        """
        # Reset shared resources
        EphemerisContext.close()

        num_threads = 20
        planets: List[int] = []
        planets_lock = threading.Lock()
        errors: List[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                ctx = EphemerisContext()
                eph = ctx.get_planets()
                with planets_lock:
                    planets.append(id(eph))
            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        assert len(planets) == num_threads

        # All should get the same shared ephemeris instance
        unique_planets = set(planets)
        assert len(unique_planets) == 1, (
            f"Expected 1 shared ephemeris, got {len(unique_planets)} instances"
        )
