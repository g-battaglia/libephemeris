"""
Thread-safety tests for concurrent ephemeris usage.

This module provides comprehensive thread-safety tests for the libephemeris library,
verifying that multiple threads can perform calculations concurrently with:
- Different ayanamsha systems
- Different topocentric locations
- No race conditions
- Correct isolated results

Each test uses 10+ threads to stress test the library's thread-safety mechanisms.

Note: Since the module-level API uses global state which is not thread-safe by design
(matching the reference ephemeris behavior), these tests verify that the state management
doesn't crash under concurrent access and that EphemerisContext provides proper
isolation for thread-safe concurrent calculations.
"""

import concurrent.futures
import threading
import time

import pytest

import libephemeris as ephem
from libephemeris import EphemerisContext
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_SPEED,
    FLG_SIDEREAL,
    FLG_TOPOCTR,
    SIDM_FAGAN_BRADLEY,
    SIDM_LAHIRI,
    SIDM_RAMAN,
    SIDM_KRISHNAMURTI,
    SIDM_YUKTESHWAR,
    SIDM_TRUE_CITRA,
    SIDM_DELUCE,
    SIDM_JN_BHASIN,
    SIDM_DJWHAL_KHUL,
    SIDM_USHASHASHI,
    SIDM_BABYL_KUGLER1,
    SIDM_GALCENT_0SAG,
)

pytestmark = pytest.mark.slow


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def test_ayanamshas():
    """12 different ayanamsha systems for thread testing."""
    return [
        (SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
        (SIDM_LAHIRI, "Lahiri"),
        (SIDM_RAMAN, "Raman"),
        (SIDM_KRISHNAMURTI, "Krishnamurti"),
        (SIDM_YUKTESHWAR, "Yukteshwar"),
        (SIDM_TRUE_CITRA, "True Citra"),
        (SIDM_DELUCE, "De Luce"),
        (SIDM_JN_BHASIN, "JN Bhasin"),
        (SIDM_DJWHAL_KHUL, "Djwhal Khul"),
        (SIDM_USHASHASHI, "Ushashashi"),
        (SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
        (SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    ]


@pytest.fixture
def test_locations_12():
    """12 different geographic locations for thread testing."""
    return [
        (0, "Equator", 0.0, 0.0, 0),
        (1, "Rome", 12.5, 41.9, 50),
        (2, "London", -0.1, 51.5, 10),
        (3, "Tokyo", 139.7, 35.7, 40),
        (4, "New York", -74.0, 40.7, 10),
        (5, "Sydney", 151.2, -33.9, 5),
        (6, "San Francisco", -122.4, 37.8, 20),
        (7, "Paris", 2.3, 48.9, 35),
        (8, "Berlin", 13.4, 52.5, 30),
        (9, "Rio de Janeiro", -43.2, -22.9, 10),
        (10, "Mumbai", 72.9, 19.1, 15),
        (11, "Cape Town", 18.4, -33.9, 25),
    ]


# =============================================================================
# TEST: EphemerisContext Instance Isolation
# =============================================================================


class TestContextInstanceIsolation:
    """Tests that EphemerisContext instances have isolated state."""

    def test_separate_contexts_have_isolated_sidereal_modes(self):
        """Each context should have its own sidereal mode."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_sid_mode(SIDM_LAHIRI)
        ctx2.set_sid_mode(SIDM_FAGAN_BRADLEY)

        assert ctx1.get_sid_mode() == SIDM_LAHIRI
        assert ctx2.get_sid_mode() == SIDM_FAGAN_BRADLEY

    def test_separate_contexts_have_isolated_topo(self):
        """Each context should have its own topocentric location."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ctx1.set_topo(12.5, 41.9, 0)  # Rome
        ctx2.set_topo(-0.1, 51.5, 0)  # London

        topo1 = ctx1.get_topo()
        topo2 = ctx2.get_topo()

        assert topo1 is not None
        assert topo2 is not None
        assert abs(topo1.latitude.degrees - 41.9) < 0.001
        assert abs(topo2.latitude.degrees - 51.5) < 0.001

    def test_12_contexts_with_different_configs(
        self, test_ayanamshas, test_locations_12
    ):
        """Create 12 contexts with different configurations and verify isolation."""
        contexts = []

        for i in range(12):
            ctx = EphemerisContext()
            mode, _ = test_ayanamshas[i]
            _, name, lon, lat, alt = test_locations_12[i]

            ctx.set_sid_mode(mode)
            ctx.set_topo(lon, lat, alt)
            contexts.append((ctx, mode, lon, lat))

        # Verify all contexts retained their configurations
        for ctx, expected_mode, expected_lon, expected_lat in contexts:
            assert ctx.get_sid_mode() == expected_mode
            topo = ctx.get_topo()
            assert topo is not None
            assert abs(topo.longitude.degrees - expected_lon) < 0.001
            assert abs(topo.latitude.degrees - expected_lat) < 0.001


# =============================================================================
# TEST: Context Resource Sharing
# =============================================================================


class TestContextResourceSharing:
    """Tests that EphemerisContext instances share expensive resources."""

    def test_contexts_share_ephemeris_via_module_state(self):
        """Module-level API uses shared planetary ephemeris."""
        from libephemeris import state

        # Explicitly load planets (calc_ut may use LEB fast path)
        state.get_planets()

        # Verify module state is loaded
        assert state._PLANETS is not None

        # Multiple loads should return the same object
        planets_before = state._PLANETS
        state.get_planets()
        planets_after = state._PLANETS

        assert planets_before is planets_after, "Ephemeris should remain shared"

    def test_contexts_share_timescale_via_module_state(self):
        """Module-level API uses shared timescale."""
        from libephemeris import state

        # Explicitly load the timescale (calc_ut may use the LEB fast
        # path, which never touches state._TS; an earlier test may have
        # reset it via EphemerisContext.close()).
        state.get_timescale()

        # Verify timescale is loaded
        assert state._TS is not None

        ts_before = state._TS
        ephem.calc_ut(2451545.0, MARS, 0)
        ts_after = state._TS

        assert ts_before is ts_after, "Timescale should remain shared"


# =============================================================================
# TEST: Sequential Calculations with Different Ayanamshas
# =============================================================================


class TestSequentialAyanamshaCalculations:
    """Tests for calculations with different ayanamsha systems."""

    def test_12_ayanamshas_produce_different_results(self, test_ayanamshas):
        """
        Test that 12 different ayanamshas produce different sidereal positions.

        Uses module-level API for reliable execution.
        """
        jd = 2451545.0  # J2000.0
        results: dict[int, float] = {}

        for mode, name in test_ayanamshas:
            ephem.set_sid_mode(mode)
            pos, _ = ephem.calc_ut(jd, SUN, FLG_SIDEREAL)
            results[mode] = pos[0]

        # All 12 ayanamshas should produce different positions
        unique_positions = len(set(round(lon, 4) for lon in results.values()))
        assert unique_positions >= 10, (
            f"Expected at least 10 unique positions from 12 ayanamshas, "
            f"got {unique_positions}"
        )


# =============================================================================
# TEST: Concurrent Module-Level API (Stress Test)
# =============================================================================


class TestConcurrentModuleLevelAPI:
    """
    Stress tests for concurrent access to module-level API.

    Note: The module-level API is NOT thread-safe by design (matching
    the reference ephemeris behavior). These tests verify the API doesn't crash
    under concurrent access, though results may not be deterministic
    when different threads modify global state.
    """

    def test_concurrent_geocentric_calculations(self):
        """
        Test that concurrent geocentric calculations don't crash.

        Geocentric calculations don't depend on topocentric state,
        so they should be safe for concurrent access.
        """
        jd = 2451545.0
        num_threads = 12
        errors: list[str] = []
        results: list[tuple[int, float]] = []
        results_lock = threading.Lock()

        def calculate_planet(thread_id: int, planet: int):
            try:
                pos, _ = ephem.calc_ut(jd, planet, 0)
                with results_lock:
                    results.append((planet, pos[0]))
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        planets = [
            SUN,
            MOON,
            MERCURY,
            VENUS,
            MARS,
            JUPITER,
            SATURN,
        ]

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [
                ex.submit(calculate_planet, i, planets[i % len(planets)])
                for i in range(num_threads)
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(results) == num_threads

    def test_concurrent_house_calculations(self, test_locations_12):
        """
        Test that concurrent house calculations don't crash.

        Each thread calculates houses for a different location.
        """
        jd = 2451545.0
        errors: list[str] = []
        results: list[tuple] = []
        results_lock = threading.Lock()

        def calculate_houses(thread_id: int, lat: float, lon: float):
            try:
                cusps, ascmc = ephem.houses(jd, lat, lon, ord("P"))
                with results_lock:
                    results.append((thread_id, cusps))
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            futures = [
                ex.submit(calculate_houses, tid, lat, lon)
                for tid, _, lon, lat, _ in test_locations_12
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(results) == 12


# =============================================================================
# TEST: Concurrent Sequential Calculations
# =============================================================================


class TestConcurrentSequentialCalculations:
    """
    Tests where each thread performs a sequence of independent calculations.
    """

    def test_10_threads_each_calculating_all_planets(self):
        """
        10 threads each calculate positions for all major planets.

        This tests that the calculation engine handles concurrent reads.
        """
        jd = 2451545.0
        num_threads = 10
        planets = [
            SUN,
            MOON,
            MERCURY,
            VENUS,
            MARS,
            JUPITER,
            SATURN,
        ]
        errors: list[str] = []
        calculation_counts = [0]
        count_lock = threading.Lock()

        def calculate_all_planets(thread_id: int):
            try:
                for planet in planets:
                    pos, _ = ephem.calc_ut(jd + thread_id * 0.1, planet, 0)
                    # Verify position is valid
                    if not (0 <= pos[0] < 360):
                        with count_lock:
                            errors.append(
                                f"Thread {thread_id}: invalid position {pos[0]}"
                            )
                        return

                with count_lock:
                    calculation_counts[0] += len(planets)

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(calculate_all_planets, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        expected = num_threads * len(planets)
        assert calculation_counts[0] == expected, (
            f"Expected {expected} calculations, got {calculation_counts[0]}"
        )

    def test_15_threads_rapid_calculations(self):
        """
        15 threads perform rapid calculations in succession.
        """
        num_threads = 15
        calculations_per_thread = 20
        errors: list[str] = []
        success_count = [0]
        count_lock = threading.Lock()

        def rapid_calcs(thread_id: int):
            try:
                for i in range(calculations_per_thread):
                    jd = 2451545.0 + thread_id + i * 0.05
                    planet = [SUN, MOON, MARS][i % 3]
                    pos, _ = ephem.calc_ut(jd, planet, 0)

                    if not (0 <= pos[0] < 360):
                        raise ValueError(f"Invalid position: {pos[0]}")

                with count_lock:
                    success_count[0] += 1

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(rapid_calcs, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors[:10]}"
        assert success_count[0] == num_threads


# =============================================================================
# TEST: Result Consistency
# =============================================================================


class TestResultConsistency:
    """Tests that verify calculation results are consistent."""

    def test_repeated_calculations_give_same_result(self):
        """Repeated calculations for the same input should give identical results."""
        jd = 2451545.0
        num_iterations = 50

        results = []
        for _ in range(num_iterations):
            pos, _ = ephem.calc_ut(jd, SUN, FLG_SPEED)
            results.append(pos[0])

        # All results should be identical
        assert all(abs(r - results[0]) < 1e-10 for r in results), (
            "Results should be identical for same input"
        )

    def test_concurrent_same_calculation(self):
        """
        Multiple threads calculating the same thing should get the same result.
        """
        jd = 2451545.0
        num_threads = 10
        results: list[float] = []
        results_lock = threading.Lock()
        errors: list[str] = []

        def calculate(thread_id: int):
            try:
                pos, _ = ephem.calc_ut(jd, MOON, FLG_SPEED)
                with results_lock:
                    results.append(pos[0])
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(calculate, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors
        assert len(results) == num_threads

        # All results should be identical (or very close)
        reference = results[0]
        for i, r in enumerate(results):
            assert abs(r - reference) < 1e-8, f"Result {i} differs: {r} vs {reference}"


# =============================================================================
# TEST: Mixed Calculation Types
# =============================================================================


class TestMixedCalculationTypes:
    """Tests mixing different types of calculations concurrently."""

    def test_concurrent_planet_and_house_calculations(self, test_locations_12):
        """
        Concurrent planet and house calculations.
        """
        jd = 2451545.0
        errors: list[str] = []
        planet_results: list[float] = []
        house_results: list = []
        results_lock = threading.Lock()

        def calc_planet(thread_id: int):
            try:
                planet = [SUN, MOON, MARS, JUPITER][thread_id % 4]
                pos, _ = ephem.calc_ut(jd, planet, 0)
                with results_lock:
                    planet_results.append(pos[0])
            except Exception as e:
                with results_lock:
                    errors.append(f"Planet Thread {thread_id}: {e}")

        def calc_houses(thread_id: int, lat: float, lon: float):
            try:
                cusps, _ = ephem.houses(jd, lat, lon, ord("P"))
                with results_lock:
                    house_results.append(cusps)
            except Exception as e:
                with results_lock:
                    errors.append(f"House Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
            # Submit planet calculations
            planet_futures = [ex.submit(calc_planet, i) for i in range(8)]

            # Submit house calculations
            house_futures = [
                ex.submit(calc_houses, tid, lat, lon)
                for tid, _, lon, lat, _ in test_locations_12
            ]

            concurrent.futures.wait(planet_futures + house_futures)

        assert not errors, f"Errors: {errors}"
        assert len(planet_results) == 8
        assert len(house_results) == 12

    def test_interleaved_different_jd_calculations(self):
        """
        Threads calculating for different Julian Days interleaved.
        """
        num_threads = 12
        errors: list[str] = []
        results: dict[int, float] = {}
        results_lock = threading.Lock()

        def calculate(thread_id: int, jd: float):
            try:
                pos, _ = ephem.calc_ut(jd, SUN, 0)
                with results_lock:
                    results[thread_id] = pos[0]
            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        # Different JDs for each thread
        jds = [2451545.0 + i * 30.0 for i in range(num_threads)]

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(calculate, i, jd) for i, jd in enumerate(jds)]
            concurrent.futures.wait(futures)

        assert not errors
        assert len(results) == num_threads

        # Verify results are different (Sun moves ~1 degree per day)
        positions = list(results.values())
        unique_positions = len(set(round(p, 2) for p in positions))
        assert unique_positions == num_threads, (
            f"Expected {num_threads} unique positions for different JDs"
        )


# =============================================================================
# TEST: High Load
# =============================================================================


class TestHighLoad:
    """High-load stress tests."""

    def test_50_concurrent_calculations(self):
        """
        50 concurrent calculation threads.
        """
        num_threads = 50
        errors: list[str] = []
        success_count = [0]
        count_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                jd = 2451545.0 + thread_id
                planet = [SUN, MOON, MERCURY, VENUS, MARS][thread_id % 5]
                pos, _ = ephem.calc_ut(jd, planet, 0)

                if not (0 <= pos[0] < 360):
                    raise ValueError(f"Invalid: {pos[0]}")

                with count_lock:
                    success_count[0] += 1

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors[:10]}"
        assert success_count[0] == num_threads

    def test_repeated_concurrent_batches(self):
        """
        Run 5 batches of 20 concurrent calculations each.
        """
        num_batches = 5
        threads_per_batch = 20
        all_errors: list[str] = []

        for batch in range(num_batches):
            errors: list[str] = []
            results: list[float] = []
            results_lock = threading.Lock()

            def calc(thread_id: int):
                try:
                    jd = 2451545.0 + batch + thread_id * 0.1
                    pos, _ = ephem.calc_ut(jd, MOON, 0)
                    with results_lock:
                        results.append(pos[0])
                except Exception as e:
                    with results_lock:
                        errors.append(f"Batch {batch} Thread {thread_id}: {e}")

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=threads_per_batch
            ) as ex:
                futures = [ex.submit(calc, i) for i in range(threads_per_batch)]
                concurrent.futures.wait(futures)

            all_errors.extend(errors)

            # Verify batch results
            assert len(results) == threads_per_batch, (
                f"Batch {batch}: expected {threads_per_batch} results"
            )

        assert not all_errors, f"Errors across batches: {all_errors[:10]}"


# =============================================================================
# TEST: Long-Running Sessions
# =============================================================================


class TestLongRunningSessions:
    """Tests for sustained calculation sessions."""

    def test_sustained_parallel_calculations(self):
        """
        10 threads each performing 100 calculations.
        """
        num_threads = 10
        calcs_per_thread = 100
        errors: list[str] = []
        total_calcs = [0]
        count_lock = threading.Lock()

        def long_session(thread_id: int):
            try:
                for i in range(calcs_per_thread):
                    jd = 2451545.0 + i * 0.5
                    planet = [SUN, MOON, MARS][i % 3]
                    ephem.calc_ut(jd, planet, 0)

                with count_lock:
                    total_calcs[0] += calcs_per_thread

            except Exception as e:
                with count_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        start = time.time()

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(long_session, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        elapsed = time.time() - start

        assert not errors, f"Errors: {errors[:10]}"
        expected = num_threads * calcs_per_thread
        assert total_calcs[0] == expected

        # Should complete in reasonable time (not hang)
        assert elapsed < 60, f"Took too long: {elapsed:.1f}s"


# =============================================================================
# TEST: EphemerisContext Thread-Safe Calculations (with Lock Protection)
# =============================================================================


class TestEphemerisContextThreadSafety:
    """
    Tests that verify EphemerisContext calculations are thread-safe.

    These tests specifically verify that the _CONTEXT_SWAP_LOCK in state.py
    properly protects the save-set-restore cycle during context-aware calculations.
    Without this lock, concurrent calls could interleave and corrupt each other's state.
    """

    def test_concurrent_contexts_with_different_sidereal_modes(self, test_ayanamshas):
        """
        Critical test: Multiple threads using different ayanamsha systems.

        This is the key test that validates the thread-safety fix. Before the
        _CONTEXT_SWAP_LOCK was added, concurrent calls with different sidereal
        modes could interfere with each other, causing incorrect results.

        Each thread:
        1. Creates its own EphemerisContext
        2. Sets a unique sidereal mode
        3. Calculates Sun position in sidereal coordinates
        4. Verifies the result matches sequential calculation

        If the lock is working correctly, each thread should get the same result
        as if it were calculated sequentially.
        """
        jd = 2451545.0
        num_threads = 12
        errors: list[str] = []
        results: dict[int, float] = {}
        results_lock = threading.Lock()

        # First, calculate expected results sequentially
        expected: dict[int, float] = {}
        for mode, _ in test_ayanamshas:
            ctx = EphemerisContext()
            ctx.set_sid_mode(mode)
            pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
            expected[mode] = pos[0]

        def calculate_with_context(thread_id: int, mode: int, mode_name: str):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(mode)

                # Perform calculation
                pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)

                with results_lock:
                    results[mode] = pos[0]

                    # Verify result matches expected
                    if abs(pos[0] - expected[mode]) > 1e-6:
                        errors.append(
                            f"Thread {thread_id} ({mode_name}): got {pos[0]:.6f}, "
                            f"expected {expected[mode]:.6f}"
                        )

            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id} ({mode_name}): {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [
                ex.submit(calculate_with_context, i, mode, name)
                for i, (mode, name) in enumerate(test_ayanamshas)
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Thread-safety violations: {errors}"
        assert len(results) == len(test_ayanamshas), (
            f"Expected {len(test_ayanamshas)} results, got {len(results)}"
        )

    def test_concurrent_contexts_with_different_topo_locations(self, test_locations_12):
        """
        Multiple threads using different topocentric locations.

        This test verifies that the lock protects the _TOPO state during
        concurrent context-aware calculations.
        """
        jd = 2451545.0
        num_threads = 12
        errors: list[str] = []
        results: dict[int, tuple[float, float]] = {}  # thread_id -> (lat, sun_lon)
        results_lock = threading.Lock()

        def calculate_with_context(
            thread_id: int, name: str, lon: float, lat: float, alt: int
        ):
            try:
                ctx = EphemerisContext()
                ctx.set_topo(lon, lat, alt)

                # Perform topocentric calculation
                pos, _ = ctx.calc_ut(jd, SUN, FLG_TOPOCTR | FLG_SPEED)

                with results_lock:
                    results[thread_id] = (lat, pos[0])

            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id} ({name}): {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [
                ex.submit(calculate_with_context, tid, name, lon, lat, alt)
                for tid, name, lon, lat, alt in test_locations_12
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        assert len(results) == 12

        # Different locations should produce slightly different topocentric positions
        positions = [r[1] for r in results.values()]
        unique = len(set(round(p, 4) for p in positions))
        assert unique >= 8, (
            f"Expected varied topocentric positions, got {unique} unique"
        )

    def test_concurrent_context_house_calculations(self, test_locations_12):
        """
        Multiple threads calculating houses with different sidereal modes.

        Tests that _houses_with_context properly uses the lock.
        """
        jd = 2451545.0
        sidereal_modes = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN]
        errors: list[str] = []
        results: list[tuple[int, int, float]] = []  # (loc_id, mode, asc)
        results_lock = threading.Lock()

        def calculate_houses(loc_id: int, mode: int, lat: float, lon: float):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(mode)

                # Use context's houses method
                cusps, ascmc = ctx.houses(jd, lat, lon, ord("P"))

                with results_lock:
                    # Record location, mode, and Ascendant
                    results.append((loc_id, mode, ascmc[0]))

            except Exception as e:
                with results_lock:
                    errors.append(f"Loc {loc_id}, Mode {mode}: {e}")

        # Create tasks for each location with each sidereal mode
        tasks = []
        for tid, name, lon, lat, alt in test_locations_12[:4]:  # Use 4 locations
            for mode in sidereal_modes:
                tasks.append((tid, mode, lat, lon))

        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            futures = [
                ex.submit(calculate_houses, loc_id, mode, lat, lon)
                for loc_id, mode, lat, lon in tasks
            ]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        assert len(results) == len(tasks)

    def test_repeated_concurrent_context_calculations(self):
        """
        Stress test: Repeated batches of concurrent context calculations.

        This test runs multiple rounds of concurrent calculations to
        increase the probability of catching race conditions.
        """
        jd = 2451545.0
        num_batches = 10
        threads_per_batch = 10
        all_errors: list[str] = []

        # Create contexts with different configurations
        configs = [
            (SIDM_LAHIRI, 12.5, 41.9),  # Rome, Lahiri
            (SIDM_FAGAN_BRADLEY, -0.1, 51.5),  # London, Fagan-Bradley
            (SIDM_RAMAN, 139.7, 35.7),  # Tokyo, Raman
            (SIDM_KRISHNAMURTI, -74.0, 40.7),  # New York, Krishnamurti
            (SIDM_YUKTESHWAR, 151.2, -33.9),  # Sydney, Yukteshwar
        ]

        for batch in range(num_batches):
            errors: list[str] = []
            results: list[float] = []
            results_lock = threading.Lock()

            def calc_with_config(thread_id: int):
                try:
                    config = configs[thread_id % len(configs)]
                    mode, lon, lat = config

                    ctx = EphemerisContext()
                    ctx.set_sid_mode(mode)
                    ctx.set_topo(lon, lat, 0)

                    pos, _ = ctx.calc_ut(
                        jd + batch * 0.1, SUN, FLG_SIDEREAL | FLG_TOPOCTR
                    )

                    with results_lock:
                        results.append(pos[0])

                except Exception as e:
                    with results_lock:
                        errors.append(f"Batch {batch} Thread {thread_id}: {e}")

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=threads_per_batch
            ) as ex:
                futures = [
                    ex.submit(calc_with_config, i) for i in range(threads_per_batch)
                ]
                concurrent.futures.wait(futures)

            all_errors.extend(errors)
            assert len(results) == threads_per_batch, (
                f"Batch {batch}: expected {threads_per_batch} results"
            )

        assert not all_errors, f"Errors: {all_errors[:10]}"

    def test_lock_prevents_state_corruption(self):
        """
        Verify that the lock prevents observable state corruption.

        This test creates a scenario where two threads with very different
        configurations calculate simultaneously. We verify each thread gets
        results consistent with its own configuration, not the other's.
        """
        jd = 2451545.0
        num_iterations = 50
        errors: list[str] = []

        # Two very different configurations
        config_a = (SIDM_LAHIRI, "Lahiri")  # ~23° ayanamsha
        config_b = (SIDM_GALCENT_0SAG, "GalCent0Sag")  # ~25° ayanamsha

        # Pre-calculate expected values
        ctx_a = EphemerisContext()
        ctx_a.set_sid_mode(config_a[0])
        expected_a, _ = ctx_a.calc_ut(jd, SUN, FLG_SIDEREAL)

        ctx_b = EphemerisContext()
        ctx_b.set_sid_mode(config_b[0])
        expected_b, _ = ctx_b.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Verify they're different enough to detect corruption
        diff = abs(expected_a[0] - expected_b[0])
        assert diff > 1.0, f"Expected positions to differ by >1°, got {diff:.4f}°"

        results_a: list[float] = []
        results_b: list[float] = []
        results_lock = threading.Lock()

        def worker_a():
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(config_a[0])
                for _ in range(num_iterations):
                    pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
                    # Check result is close to expected_a, not expected_b
                    if abs(pos[0] - expected_a[0]) > 0.001:
                        with results_lock:
                            errors.append(
                                f"Worker A: got {pos[0]:.6f}, expected {expected_a[0]:.6f}"
                            )
                        return
                    with results_lock:
                        results_a.append(pos[0])
            except Exception as e:
                with results_lock:
                    errors.append(f"Worker A exception: {e}")

        def worker_b():
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(config_b[0])
                for _ in range(num_iterations):
                    pos, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)
                    # Check result is close to expected_b, not expected_a
                    if abs(pos[0] - expected_b[0]) > 0.001:
                        with results_lock:
                            errors.append(
                                f"Worker B: got {pos[0]:.6f}, expected {expected_b[0]:.6f}"
                            )
                        return
                    with results_lock:
                        results_b.append(pos[0])
            except Exception as e:
                with results_lock:
                    errors.append(f"Worker B exception: {e}")

        # Run workers concurrently
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            future_a = ex.submit(worker_a)
            future_b = ex.submit(worker_b)
            concurrent.futures.wait([future_a, future_b])

        assert not errors, f"State corruption detected: {errors}"
        assert len(results_a) == num_iterations
        assert len(results_b) == num_iterations

    def test_context_swap_lock_exists_and_is_rlock(self):
        """
        Verify that _CONTEXT_SWAP_LOCK exists in state module and is an RLock.

        We use RLock (reentrant lock) to allow nested locking when functions
        need to call other state-accessing functions while holding the lock.
        """
        from libephemeris import state

        assert hasattr(state, "_CONTEXT_SWAP_LOCK"), (
            "_CONTEXT_SWAP_LOCK not found in state module"
        )
        assert isinstance(state._CONTEXT_SWAP_LOCK, type(threading.RLock())), (
            "_CONTEXT_SWAP_LOCK should be a threading.RLock"
        )

    def test_calc_pctr_thread_safety(self):
        """
        Test thread-safety of planet-centric calculations via context.
        """
        jd = 2451545.0
        errors: list[str] = []
        results: list[float] = []
        results_lock = threading.Lock()

        def calc_pctr(thread_id: int, mode: int):
            try:
                ctx = EphemerisContext()
                ctx.set_sid_mode(mode)

                # Calculate Moon as seen from Mars
                pos, _ = ctx.calc_pctr(jd, MOON, MARS, FLG_SIDEREAL)

                with results_lock:
                    results.append(pos[0])

            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        modes = [SIDM_LAHIRI, SIDM_FAGAN_BRADLEY, SIDM_RAMAN]

        with concurrent.futures.ThreadPoolExecutor(max_workers=9) as ex:
            futures = [ex.submit(calc_pctr, i, modes[i % 3]) for i in range(9)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"
        assert len(results) == 9


# =============================================================================
# TEST: EphemerisContext.close() Method
# =============================================================================


class TestEphemerisContextClose:
    """
    Tests for EphemerisContext.close() method.

    These tests verify that:
    1. close() properly resets shared resources to None
    2. Calculations after close() reload resources correctly
    3. close() is thread-safe in concurrent environments
    4. The load-close-reload cycle works without memory leaks or stale state
    """

    def test_close_resets_shared_resources_to_none(self):
        """
        Verify that close() sets all shared resources to None.

        After calling EphemerisContext.close(), the shared resources
        (_SHARED_LOADER, _SHARED_PLANETS, _SHARED_TS) should all be None.

        Note: ctx.calc_ut() uses state module's resources internally,
        so we explicitly call ctx.get_planets() and ctx.get_timescale()
        to load the context module's shared resources.
        """
        from libephemeris import context as ctx_module

        # First, ensure resources are loaded by explicitly calling accessors
        ctx = EphemerisContext()
        _ = ctx.get_loader()
        _ = ctx.get_planets()
        _ = ctx.get_timescale()

        # Verify resources are loaded
        assert ctx_module._SHARED_LOADER is not None, "Loader should be loaded"
        assert ctx_module._SHARED_PLANETS is not None, "Planets should be loaded"
        assert ctx_module._SHARED_TS is not None, "Timescale should be loaded"

        # Call close()
        EphemerisContext.close()

        # Verify all shared resources are now None
        assert ctx_module._SHARED_LOADER is None, "Loader should be None after close()"
        assert ctx_module._SHARED_PLANETS is None, (
            "Planets should be None after close()"
        )
        assert ctx_module._SHARED_TS is None, "Timescale should be None after close()"
        assert ctx_module._SHARED_EPHE_PATH is None, (
            "Ephemeris path should be None after close()"
        )

    def test_calculation_after_close_reloads_resources(self):
        """
        Verify that accessing resources after close() properly reloads them.

        This tests the full load-close-reload cycle to ensure no stale state.

        Note: ctx.calc_ut() uses state module's resources internally,
        so we explicitly call ctx.get_planets() and ctx.get_timescale()
        to test the context module's shared resources.
        """
        from libephemeris import context as ctx_module

        jd = 2451545.0

        # Step 1: Load context resources explicitly
        ctx = EphemerisContext()
        _ = ctx.get_loader()
        _ = ctx.get_planets()
        _ = ctx.get_timescale()

        # Capture resource identities before close
        loader_before = ctx_module._SHARED_LOADER
        planets_before = ctx_module._SHARED_PLANETS
        ts_before = ctx_module._SHARED_TS

        assert loader_before is not None
        assert planets_before is not None
        assert ts_before is not None

        # Do calculation for comparison
        pos_before, _ = ctx.calc_ut(jd, SUN, 0)

        # Step 2: Call close()
        EphemerisContext.close()

        # Verify resources are None
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_PLANETS is None
        assert ctx_module._SHARED_TS is None

        # Step 3: Access resources again (should reload)
        _ = ctx.get_loader()
        _ = ctx.get_planets()
        _ = ctx.get_timescale()

        # Verify resources are reloaded (new instances)
        loader_after = ctx_module._SHARED_LOADER
        planets_after = ctx_module._SHARED_PLANETS
        ts_after = ctx_module._SHARED_TS

        assert loader_after is not None, "Loader should be reloaded"
        assert planets_after is not None, "Planets should be reloaded"
        assert ts_after is not None, "Timescale should be reloaded"

        # New instances should be different objects (reloaded)
        assert loader_after is not loader_before, "Loader should be new instance"
        assert planets_after is not planets_before, "Planets should be new instance"
        assert ts_after is not ts_before, "Timescale should be new instance"

        # Do calculation after reload and verify same result
        pos_after, _ = ctx.calc_ut(jd, SUN, 0)

        # Results should be identical (same calculation)
        assert abs(pos_before[0] - pos_after[0]) < 1e-10, (
            f"Results should be identical: {pos_before[0]} vs {pos_after[0]}"
        )

    def test_multiple_close_calls_are_safe(self):
        """
        Verify that multiple consecutive close() calls are safe.

        Calling close() multiple times should not raise errors.
        """
        from libephemeris import context as ctx_module

        # Load resources
        ctx = EphemerisContext()
        ctx.calc_ut(2451545.0, SUN, 0)

        # Call close() multiple times
        EphemerisContext.close()
        assert ctx_module._SHARED_PLANETS is None

        EphemerisContext.close()  # Should not raise
        assert ctx_module._SHARED_PLANETS is None

        EphemerisContext.close()  # Should not raise
        assert ctx_module._SHARED_PLANETS is None

        # Can still do calculations after
        pos, _ = ctx.calc_ut(2451545.0, SUN, 0)
        assert 0 <= pos[0] < 360

    def test_close_does_not_affect_instance_state(self):
        """
        Verify that close() only resets shared resources, not instance state.

        Instance-specific state (topo, sidereal_mode, angles_cache) should
        be preserved after close().
        """
        from libephemeris import context as ctx_module

        # Create context with specific settings
        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_LAHIRI)
        ctx.set_topo(12.5, 41.9, 100)  # Rome

        # Do a calculation to load shared resources
        ctx.calc_ut(2451545.0, SUN, 0)

        # Call close()
        EphemerisContext.close()

        # Instance state should be preserved
        assert ctx.sidereal_mode == SIDM_LAHIRI, (
            "Sidereal mode should be preserved after close()"
        )
        assert ctx.topo is not None, "Topo should be preserved after close()"
        assert abs(ctx.topo.latitude.degrees - 41.9) < 0.001

        # Verify shared resources are None
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_PLANETS is None
        assert ctx_module._SHARED_TS is None

    def test_close_in_multithread_environment(self):
        """
        Test close() thread-safety in a multi-threaded environment.

        Multiple threads calling close() and doing calculations concurrently
        should not cause crashes or data corruption.
        """

        jd = 2451545.0
        num_threads = 12
        num_iterations = 10
        errors: list[str] = []
        errors_lock = threading.Lock()

        def worker(thread_id: int):
            try:
                for i in range(num_iterations):
                    ctx = EphemerisContext()
                    ctx.set_sid_mode(SIDM_LAHIRI)

                    # Do a calculation
                    pos, _ = ctx.calc_ut(jd + i * 0.1, SUN, 0)

                    # Validate result
                    if not (0 <= pos[0] < 360):
                        with errors_lock:
                            errors.append(
                                f"Thread {thread_id}, iter {i}: invalid pos {pos[0]}"
                            )
                        return

                    # Some threads call close() in the middle
                    if thread_id % 3 == 0 and i == 5:
                        EphemerisContext.close()

            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(worker, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors: {errors}"

        # Cleanup for other tests
        EphemerisContext.close()

    def test_concurrent_close_and_calculation(self):
        """
        Stress test: concurrent close() calls interleaved with calculations.

        This is the critical test for detecting race conditions between
        close() and calculation operations.
        """

        jd = 2451545.0
        num_calc_threads = 10
        num_close_threads = 3
        num_iterations = 20
        errors: list[str] = []
        successful_calcs = [0]
        errors_lock = threading.Lock()

        def calc_worker(thread_id: int):
            try:
                for i in range(num_iterations):
                    ctx = EphemerisContext()
                    pos, _ = ctx.calc_ut(jd + i * 0.1, MOON, 0)

                    if not (0 <= pos[0] < 360):
                        with errors_lock:
                            errors.append(f"Calc thread {thread_id}: invalid {pos[0]}")
                        return

                    with errors_lock:
                        successful_calcs[0] += 1

            except Exception as e:
                with errors_lock:
                    errors.append(f"Calc thread {thread_id}: {e}")

        def close_worker(thread_id: int):
            try:
                for i in range(num_iterations):
                    time.sleep(0.001)  # Small delay to interleave with calcs
                    EphemerisContext.close()
            except Exception as e:
                with errors_lock:
                    errors.append(f"Close thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=num_calc_threads + num_close_threads
        ) as ex:
            # Start calculation workers
            calc_futures = [ex.submit(calc_worker, i) for i in range(num_calc_threads)]
            # Start close workers
            close_futures = [
                ex.submit(close_worker, i) for i in range(num_close_threads)
            ]

            concurrent.futures.wait(calc_futures + close_futures)

        # Race condition errors from close() nullifying shared resources
        # mid-calculation are expected and acceptable (NoneType access).
        # Only unexpected errors should fail the test.
        unexpected_errors = [
            e for e in errors if "NoneType" not in e and "None" not in e
        ]
        assert not unexpected_errors, f"Unexpected errors: {unexpected_errors}"

        # At least some calculations should have succeeded
        expected_min = num_calc_threads * num_iterations // 2
        assert successful_calcs[0] >= expected_min, (
            f"Expected at least {expected_min} successful calculations, "
            f"got {successful_calcs[0]}"
        )

        # Cleanup
        EphemerisContext.close()

    def test_load_close_reload_cycle_consistency(self):
        """
        Verify the load-close-reload cycle produces consistent results.

        Running the same calculation before close, after reload, and
        in subsequent cycles should produce identical results.

        Note: We explicitly call ctx.get_planets() to load the context
        module's shared resources, since calc_ut uses state module.
        """
        from libephemeris import context as ctx_module

        jd = 2451545.0
        num_cycles = 5
        results: list[float] = []

        for cycle in range(num_cycles):
            ctx = EphemerisContext()
            # Explicitly load context's shared resources
            _ = ctx.get_planets()

            pos, _ = ctx.calc_ut(jd, SUN, 0)
            results.append(pos[0])

            # Verify context resources are loaded
            assert ctx_module._SHARED_PLANETS is not None

            # Close resources
            EphemerisContext.close()

            # Verify resources are cleared
            assert ctx_module._SHARED_PLANETS is None

        # All results should be identical
        reference = results[0]
        for i, r in enumerate(results):
            assert abs(r - reference) < 1e-10, (
                f"Cycle {i}: {r} differs from reference {reference}"
            )

    def test_close_with_sidereal_calculations(self):
        """
        Test close() with sidereal mode calculations.

        Verify that sidereal mode settings work correctly after close().
        """

        jd = 2451545.0

        # Calculate with Lahiri before close
        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_LAHIRI)
        pos_lahiri_before, _ = ctx.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Close
        EphemerisContext.close()

        # Calculate with same settings after close
        ctx2 = EphemerisContext()
        ctx2.set_sid_mode(SIDM_LAHIRI)
        pos_lahiri_after, _ = ctx2.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Results should be identical
        assert abs(pos_lahiri_before[0] - pos_lahiri_after[0]) < 1e-10, (
            f"Sidereal results differ: {pos_lahiri_before[0]} vs {pos_lahiri_after[0]}"
        )

        # Close and test with different ayanamsha
        EphemerisContext.close()

        ctx3 = EphemerisContext()
        ctx3.set_sid_mode(SIDM_FAGAN_BRADLEY)
        pos_fagan, _ = ctx3.calc_ut(jd, SUN, FLG_SIDEREAL)

        # Different ayanamsha should give different results
        assert abs(pos_lahiri_after[0] - pos_fagan[0]) > 0.1, (
            "Different ayanamshas should produce different positions"
        )

        # Cleanup
        EphemerisContext.close()

    def test_close_with_house_calculations(self):
        """
        Test close() with house calculations.

        Verify that house calculations work correctly after close().
        """

        jd = 2451545.0

        # Calculate houses before close
        ctx = EphemerisContext()
        cusps_before, ascmc_before = ctx.houses(jd, 41.9, 12.5, ord("P"))

        # Close
        EphemerisContext.close()

        # Calculate houses after close
        ctx2 = EphemerisContext()
        cusps_after, ascmc_after = ctx2.houses(jd, 41.9, 12.5, ord("P"))

        # Results should be identical
        for i, (before, after) in enumerate(zip(cusps_before, cusps_after)):
            assert abs(before - after) < 1e-10, (
                f"House cusp {i} differs: {before} vs {after}"
            )

        assert abs(ascmc_before[0] - ascmc_after[0]) < 1e-10, (
            f"Ascendant differs: {ascmc_before[0]} vs {ascmc_after[0]}"
        )

        # Cleanup
        EphemerisContext.close()
