"""
Extended tests for EphemerisContext lazy initialization.

This module tests the get_loader(), get_timescale(), and get_planets() methods
which use double-checked locking for thread-safe lazy initialization of shared
resources.

These tests explicitly verify:
- Lazy initialization (calling methods before any calculation)
- Singleton behavior (multiple calls return same object)
- Resource reset after EphemerisContext.close()
- Thread-safe concurrent initialization
"""

import concurrent.futures
import threading

import pytest
from skyfield.api import Loader
from skyfield.jpllib import SpiceKernel
from skyfield.timelib import Timescale

from libephemeris import EphemerisContext
from libephemeris import context as ctx_module
from libephemeris.constants import SUN


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(autouse=True)
def reset_shared_resources():
    """Reset shared resources before and after each test."""
    # Close any existing shared resources before test
    EphemerisContext.close()
    yield
    # Clean up after test
    EphemerisContext.close()


# =============================================================================
# TEST: Lazy Initialization of get_loader()
# =============================================================================


class TestGetLoaderLazyInitialization:
    """Tests for get_loader() lazy initialization and singleton behavior."""

    def test_get_loader_before_any_calculation(self):
        """
        Verify get_loader() can be called before any calculation.

        This tests the lazy initialization path where _SHARED_LOADER is None
        and must be created.
        """
        # Ensure we start with clean state
        assert ctx_module._SHARED_LOADER is None, "Expected _SHARED_LOADER to be None"

        # Create context and call get_loader() directly
        ctx = EphemerisContext()
        loader = ctx.get_loader()

        # Verify loader is valid
        assert loader is not None
        assert isinstance(loader, Loader)

        # Verify shared state was updated
        assert ctx_module._SHARED_LOADER is not None
        assert ctx_module._SHARED_LOADER is loader

    def test_multiple_get_loader_calls_return_same_object(self):
        """
        Multiple calls to get_loader() should return the same instance.

        This verifies the singleton pattern is working correctly.
        """
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        loader1 = ctx1.get_loader()
        loader2 = ctx1.get_loader()  # Same context
        loader3 = ctx2.get_loader()  # Different context

        # All should be the exact same object
        assert loader1 is loader2, "Same context should return same loader"
        assert loader1 is loader3, "Different contexts should share the same loader"

    def test_get_loader_multiple_contexts_same_loader(self):
        """All contexts should share the same loader instance."""
        contexts = [EphemerisContext() for _ in range(10)]
        loaders = [ctx.get_loader() for ctx in contexts]

        first_loader = loaders[0]
        for i, loader in enumerate(loaders):
            assert loader is first_loader, (
                f"Context {i} returned different loader instance"
            )

    def test_get_loader_after_close_creates_new_instance(self):
        """
        After close(), get_loader() should create a new loader instance.

        This tests the resource reset behavior.
        """
        ctx = EphemerisContext()

        # Get initial loader
        loader_before = ctx.get_loader()
        assert loader_before is not None

        # Close shared resources
        EphemerisContext.close()

        # Verify shared state is reset
        assert ctx_module._SHARED_LOADER is None

        # Get loader again - should create new instance
        loader_after = ctx.get_loader()
        assert loader_after is not None
        assert isinstance(loader_after, Loader)

        # Note: We don't compare identity since a new loader with the same
        # directory may or may not be equal depending on implementation


# =============================================================================
# TEST: Lazy Initialization of get_timescale()
# =============================================================================


class TestGetTimescaleLazyInitialization:
    """Tests for get_timescale() lazy initialization and singleton behavior."""

    def test_get_timescale_before_any_calculation(self):
        """
        Verify get_timescale() can be called before any calculation.

        This tests lazy initialization where _SHARED_TS is None.
        """
        assert ctx_module._SHARED_TS is None, "Expected _SHARED_TS to be None"

        ctx = EphemerisContext()
        ts = ctx.get_timescale()

        assert ts is not None
        assert isinstance(ts, Timescale)
        assert ctx_module._SHARED_TS is not None
        assert ctx_module._SHARED_TS is ts

    def test_get_timescale_does_not_initialize_loader(self):
        """
        get_timescale() reuses state.get_timescale() (one timescale per
        process) and must NOT initialize the context's own loader as a
        side effect.
        """
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_TS is None

        ctx = EphemerisContext()
        ts = ctx.get_timescale()

        # The shared timescale is set; the loader stays untouched
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_TS is not None
        assert ts is ctx_module._SHARED_TS

    def test_multiple_get_timescale_calls_return_same_object(self):
        """Multiple calls should return the same timescale instance."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        ts1 = ctx1.get_timescale()
        ts2 = ctx1.get_timescale()
        ts3 = ctx2.get_timescale()

        assert ts1 is ts2, "Same context should return same timescale"
        assert ts1 is ts3, "Different contexts should share the same timescale"

    def test_get_timescale_after_close_creates_new_instance(self):
        """After close(), get_timescale() should create a new instance."""
        ctx = EphemerisContext()

        ts_before = ctx.get_timescale()
        assert ts_before is not None

        EphemerisContext.close()

        assert ctx_module._SHARED_TS is None
        assert ctx_module._SHARED_LOADER is None

        ts_after = ctx.get_timescale()
        assert ts_after is not None
        assert isinstance(ts_after, Timescale)


# =============================================================================
# TEST: Lazy Initialization of get_planets()
# =============================================================================


class TestGetPlanetsLazyInitialization:
    """Tests for get_planets() lazy initialization and singleton behavior."""

    def test_get_planets_before_any_calculation(self):
        """
        Verify get_planets() can be called before any calculation.

        This tests lazy initialization where _SHARED_PLANETS is None.
        """
        assert ctx_module._SHARED_PLANETS is None, "Expected _SHARED_PLANETS to be None"

        ctx = EphemerisContext()
        planets = ctx.get_planets()

        assert planets is not None
        assert isinstance(planets, SpiceKernel)
        assert ctx_module._SHARED_PLANETS is not None
        assert ctx_module._SHARED_PLANETS is planets

    def test_get_planets_also_initializes_loader(self):
        """
        get_planets() should also initialize the loader as a dependency.
        """
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_PLANETS is None

        ctx = EphemerisContext()
        planets = ctx.get_planets()

        # Both should now be initialized
        assert ctx_module._SHARED_LOADER is not None
        assert ctx_module._SHARED_PLANETS is not None
        assert planets is ctx_module._SHARED_PLANETS

    def test_multiple_get_planets_calls_return_same_object(self):
        """Multiple calls should return the same SpiceKernel instance."""
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        planets1 = ctx1.get_planets()
        planets2 = ctx1.get_planets()
        planets3 = ctx2.get_planets()

        assert planets1 is planets2, "Same context should return same planets"
        assert planets1 is planets3, "Different contexts should share the same planets"

    def test_get_planets_after_close_creates_new_instance(self):
        """After close(), get_planets() should create a new instance."""
        ctx = EphemerisContext()

        planets_before = ctx.get_planets()
        assert planets_before is not None

        EphemerisContext.close()

        assert ctx_module._SHARED_PLANETS is None
        assert ctx_module._SHARED_LOADER is None

        planets_after = ctx.get_planets()
        assert planets_after is not None
        assert isinstance(planets_after, SpiceKernel)


# =============================================================================
# TEST: Thread-Safety of Lazy Initialization
# =============================================================================


class TestLazyInitializationThreadSafety:
    """Tests for thread-safe lazy initialization of shared resources."""

    def test_concurrent_get_loader_all_return_same_instance(self):
        """
        Multiple threads calling get_loader() concurrently should all
        get the same loader instance.

        This tests the double-checked locking pattern.
        """
        num_threads = 20
        loaders: list[Loader] = []
        loaders_lock = threading.Lock()
        errors: list[str] = []

        def get_loader_in_thread(thread_id: int):
            try:
                ctx = EphemerisContext()
                loader = ctx.get_loader()
                with loaders_lock:
                    loaders.append(loader)
            except Exception as e:
                with loaders_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(get_loader_in_thread, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(loaders) == num_threads

        # All loaders should be the same instance
        first_loader = loaders[0]
        for i, loader in enumerate(loaders):
            assert loader is first_loader, f"Thread {i} got different loader instance"

    def test_concurrent_get_timescale_all_return_same_instance(self):
        """
        Multiple threads calling get_timescale() concurrently should all
        get the same timescale instance.
        """
        num_threads = 20
        timescales: list[Timescale] = []
        ts_lock = threading.Lock()
        errors: list[str] = []

        def get_ts_in_thread(thread_id: int):
            try:
                ctx = EphemerisContext()
                ts = ctx.get_timescale()
                with ts_lock:
                    timescales.append(ts)
            except Exception as e:
                with ts_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(get_ts_in_thread, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(timescales) == num_threads

        first_ts = timescales[0]
        for i, ts in enumerate(timescales):
            assert ts is first_ts, f"Thread {i} got different timescale instance"

    def test_concurrent_get_planets_all_return_same_instance(self):
        """
        Multiple threads calling get_planets() concurrently should all
        get the same SpiceKernel instance.
        """
        num_threads = 20
        kernels: list[SpiceKernel] = []
        kernels_lock = threading.Lock()
        errors: list[str] = []

        def get_planets_in_thread(thread_id: int):
            try:
                ctx = EphemerisContext()
                planets = ctx.get_planets()
                with kernels_lock:
                    kernels.append(planets)
            except Exception as e:
                with kernels_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(get_planets_in_thread, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(kernels) == num_threads

        first_kernel = kernels[0]
        for i, kernel in enumerate(kernels):
            assert kernel is first_kernel, (
                f"Thread {i} got different SpiceKernel instance"
            )

    def test_concurrent_mixed_resource_access(self):
        """
        Multiple threads accessing different resources concurrently should
        all get the shared instances without race conditions.
        """
        num_threads = 30
        results: list[tuple[str, object]] = []
        results_lock = threading.Lock()
        errors: list[str] = []

        def access_resource(thread_id: int):
            try:
                ctx = EphemerisContext()
                resource_type = thread_id % 3
                if resource_type == 0:
                    resource = ctx.get_loader()
                    name = "loader"
                elif resource_type == 1:
                    resource = ctx.get_timescale()
                    name = "timescale"
                else:
                    resource = ctx.get_planets()
                    name = "planets"

                with results_lock:
                    results.append((name, resource))

            except Exception as e:
                with results_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
            futures = [ex.submit(access_resource, i) for i in range(num_threads)]
            concurrent.futures.wait(futures)

        assert not errors, f"Errors occurred: {errors}"
        assert len(results) == num_threads

        # Group results by type and verify all same
        loaders = [r for name, r in results if name == "loader"]
        timescales = [r for name, r in results if name == "timescale"]
        kernels = [r for name, r in results if name == "planets"]

        # All loaders should be the same instance
        if loaders:
            first_loader = loaders[0]
            assert all(loader is first_loader for loader in loaders)

        # All timescales should be the same instance
        if timescales:
            first_ts = timescales[0]
            assert all(ts is first_ts for ts in timescales)

        # All kernels should be the same instance
        if kernels:
            first_kernel = kernels[0]
            assert all(k is first_kernel for k in kernels)


# =============================================================================
# TEST: Resource Lifecycle After close()
# =============================================================================


class TestResourceLifecycleAfterClose:
    """Tests for resource lifecycle management with close() method."""

    def test_close_resets_all_shared_resources(self):
        """
        close() should reset all shared resources to None.
        """
        ctx = EphemerisContext()

        # Initialize all shared resources
        _ = ctx.get_loader()
        _ = ctx.get_timescale()
        _ = ctx.get_planets()

        # Verify all are initialized
        assert ctx_module._SHARED_LOADER is not None
        assert ctx_module._SHARED_TS is not None
        assert ctx_module._SHARED_PLANETS is not None

        # Close
        EphemerisContext.close()

        # Verify all are reset
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_TS is None
        assert ctx_module._SHARED_PLANETS is None

    def test_close_resets_shared_ephemeris_config(self):
        """
        close() should reset shared ephemeris path and file to defaults.
        """
        ctx = EphemerisContext(ephe_path="/custom/path", ephe_file="custom.bsp")

        # Verify custom config was set
        assert ctx_module._SHARED_EPHE_PATH == "/custom/path"
        assert ctx_module._SHARED_EPHE_FILE == "custom.bsp"

        EphemerisContext.close()

        # Verify reset to defaults
        assert ctx_module._SHARED_EPHE_PATH is None
        assert ctx_module._SHARED_EPHE_FILE == "de440.bsp"

    def test_calculation_after_close_works(self):
        """
        After close(), calculations should still work by reinitializing resources.
        """
        ctx = EphemerisContext()

        # Do a calculation
        pos1, _ = ctx.calc_ut(2451545.0, SUN, 0)
        assert 0 <= pos1[0] < 360

        # Close
        EphemerisContext.close()

        # Verify resources are closed
        assert ctx_module._SHARED_PLANETS is None
        assert ctx_module._SHARED_TS is None

        # Do another calculation - should reinitialize automatically
        pos2, _ = ctx.calc_ut(2451545.0, SUN, 0)
        assert 0 <= pos2[0] < 360

        # Results should be consistent
        assert abs(pos1[0] - pos2[0]) < 1e-6, (
            f"Position mismatch after close: {pos1[0]} vs {pos2[0]}"
        )

    def test_close_is_idempotent(self):
        """
        Calling close() multiple times should not cause errors.
        """
        ctx = EphemerisContext()
        _ = ctx.get_planets()

        # Close multiple times - should not raise
        EphemerisContext.close()
        EphemerisContext.close()
        EphemerisContext.close()

        # Verify state is clean
        assert ctx_module._SHARED_LOADER is None
        assert ctx_module._SHARED_PLANETS is None

    def test_close_with_no_initialization(self):
        """
        Calling close() when nothing was initialized should not raise.
        """
        # Nothing is initialized (fixture runs reset)
        assert ctx_module._SHARED_LOADER is None

        # Should not raise
        EphemerisContext.close()

        assert ctx_module._SHARED_LOADER is None


# =============================================================================
# TEST: Instance State Persistence After close()
# =============================================================================


class TestInstanceStateAfterClose:
    """Tests that instance-specific state survives close()."""

    def test_topo_persists_after_close(self):
        """
        Instance topo setting should persist after close() since it's
        instance-specific, not shared.
        """
        ctx = EphemerisContext()
        ctx.set_topo(12.5, 41.9, 0)  # Rome

        topo_before = ctx.get_topo()
        assert topo_before is not None

        EphemerisContext.close()

        topo_after = ctx.get_topo()
        assert topo_after is not None
        assert topo_after is topo_before  # Same object

    def test_sidereal_mode_persists_after_close(self):
        """
        Instance sidereal_mode should persist after close().
        """
        from libephemeris.constants import SIDM_FAGAN_BRADLEY

        ctx = EphemerisContext()
        ctx.set_sid_mode(SIDM_FAGAN_BRADLEY)

        mode_before = ctx.get_sid_mode()

        EphemerisContext.close()

        mode_after = ctx.get_sid_mode()
        assert mode_after == mode_before == SIDM_FAGAN_BRADLEY

    def test_angles_cache_persists_after_close(self):
        """
        Instance angles cache should persist after close().
        """
        ctx = EphemerisContext()
        ctx.set_angles_cache({"Sun": 120.5, "Moon": 240.3})

        cache_before = ctx.get_angles_cache()

        EphemerisContext.close()

        cache_after = ctx.get_angles_cache()
        assert cache_after == cache_before
        assert cache_after["Sun"] == 120.5


# =============================================================================
# TEST: Angles Cache Isolation (set/get/clear)
# =============================================================================


class TestAnglesCacheIsolation:
    """
    Tests for angles cache isolation behavior.

    These tests verify that set_angles_cache(), get_angles_cache(), and
    clear_angles_cache() work correctly and that the cache is protected
    from external mutations via copy-on-set.

    Related to context.py:249-268.
    """

    def test_set_angles_cache_creates_copy(self):
        """
        Verify set_angles_cache() creates a copy of the input dict.

        After setting the cache, modifying the original dict should NOT
        affect the cached values. This validates the copy() at context.py:260.
        """
        ctx = EphemerisContext()

        # Create original dict
        original_angles = {"Sun": 120.5, "Moon": 240.3, "Asc": 15.7}

        # Set the cache
        ctx.set_angles_cache(original_angles)

        # Modify the original dict
        original_angles["Sun"] = 999.0
        original_angles["Mars"] = 45.0
        del original_angles["Moon"]

        # Verify cache is unaffected
        cached = ctx.get_angles_cache()
        assert cached["Sun"] == 120.5, (
            "Cache should not be affected by original dict mutation"
        )
        assert cached["Moon"] == 240.3, "Deleted key should still exist in cache"
        assert cached["Asc"] == 15.7
        assert "Mars" not in cached, "New key should not appear in cache"

    def test_get_angles_cache_after_set(self):
        """
        Verify get_angles_cache() returns the values set by set_angles_cache().

        This is a basic integration test ensuring set/get work together.
        """
        ctx = EphemerisContext()

        test_angles = {
            "Sun": 0.0,
            "Moon": 90.0,
            "Mercury": 180.0,
            "Venus": 270.0,
            "Asc": 45.5,
            "MC": 135.5,
        }

        ctx.set_angles_cache(test_angles)
        cached = ctx.get_angles_cache()

        # Verify all values are present and correct
        assert len(cached) == 6
        assert cached["Sun"] == 0.0
        assert cached["Moon"] == 90.0
        assert cached["Mercury"] == 180.0
        assert cached["Venus"] == 270.0
        assert cached["Asc"] == 45.5
        assert cached["MC"] == 135.5

    def test_clear_angles_cache_resets_to_empty_dict(self):
        """
        Verify clear_angles_cache() resets the cache to an empty dict.

        This tests the reset behavior at context.py:268.
        """
        ctx = EphemerisContext()

        # Set some values
        ctx.set_angles_cache({"Sun": 120.5, "Moon": 240.3})

        # Verify cache is populated
        assert len(ctx.get_angles_cache()) == 2

        # Clear the cache
        ctx.clear_angles_cache()

        # Verify cache is now empty
        cached = ctx.get_angles_cache()
        assert cached == {}
        assert len(cached) == 0
        assert isinstance(cached, dict)

    def test_clear_angles_cache_on_empty_cache(self):
        """
        Verify clear_angles_cache() works even when cache is already empty.

        Should not raise any errors.
        """
        ctx = EphemerisContext()

        # Cache starts empty
        assert ctx.get_angles_cache() == {}

        # Clear should not raise
        ctx.clear_angles_cache()

        # Should still be empty
        assert ctx.get_angles_cache() == {}

    def test_set_angles_cache_overwrites_previous(self):
        """
        Verify set_angles_cache() completely replaces previous cache.

        Setting a new cache should replace all previous values, not merge.
        """
        ctx = EphemerisContext()

        # Set initial cache
        ctx.set_angles_cache({"Sun": 100.0, "Moon": 200.0, "Mars": 300.0})

        # Set new cache with different keys
        ctx.set_angles_cache({"Venus": 50.0, "Jupiter": 150.0})

        cached = ctx.get_angles_cache()

        # Should only have new keys
        assert "Sun" not in cached
        assert "Moon" not in cached
        assert "Mars" not in cached
        assert cached["Venus"] == 50.0
        assert cached["Jupiter"] == 150.0
        assert len(cached) == 2

    def test_set_angles_cache_with_empty_dict(self):
        """
        Verify set_angles_cache({}) clears the cache.

        Setting an empty dict should result in an empty cache.
        """
        ctx = EphemerisContext()

        # Set initial values
        ctx.set_angles_cache({"Sun": 120.5})
        assert len(ctx.get_angles_cache()) == 1

        # Set empty dict
        ctx.set_angles_cache({})

        # Cache should be empty
        assert ctx.get_angles_cache() == {}

    def test_get_angles_cache_returns_same_internal_dict(self):
        """
        Verify multiple get_angles_cache() calls return the same internal dict.

        This tests that get_angles_cache returns the actual internal reference.
        Note: While this allows mutation, set_angles_cache creates a copy on input.
        """
        ctx = EphemerisContext()

        ctx.set_angles_cache({"Sun": 120.5})

        cache1 = ctx.get_angles_cache()
        cache2 = ctx.get_angles_cache()

        # Should be the same object
        assert cache1 is cache2

    def test_angles_cache_instance_isolation(self):
        """
        Verify each EphemerisContext instance has its own angles cache.

        Changes to one context's cache should not affect another's.
        """
        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        # Set different values in each context
        ctx1.set_angles_cache({"Sun": 100.0})
        ctx2.set_angles_cache({"Moon": 200.0})

        # Verify isolation
        cache1 = ctx1.get_angles_cache()
        cache2 = ctx2.get_angles_cache()

        assert cache1 == {"Sun": 100.0}
        assert cache2 == {"Moon": 200.0}
        assert "Moon" not in cache1
        assert "Sun" not in cache2

    def test_angles_cache_with_float_precision(self):
        """
        Verify angles cache preserves float precision.

        Arabic parts calculations require precise angles.
        """
        ctx = EphemerisContext()

        precise_angles = {
            "Sun": 123.456789012345,
            "Moon": 0.000000000001,
            "Asc": 359.999999999999,
        }

        ctx.set_angles_cache(precise_angles)
        cached = ctx.get_angles_cache()

        assert cached["Sun"] == 123.456789012345
        assert cached["Moon"] == 0.000000000001
        assert cached["Asc"] == 359.999999999999


# =============================================================================
# TEST: Concurrent close() and Initialization
# =============================================================================


class TestConcurrentCloseAndInit:
    """Tests for thread-safety when close() and initialization happen concurrently."""

    def test_close_during_concurrent_initialization(self):
        """
        If close() is called while other threads are initializing,
        the system should remain consistent (no crashes, no corruption).
        """
        num_threads = 20
        errors: list[str] = []
        errors_lock = threading.Lock()

        def init_resources(thread_id: int):
            try:
                ctx = EphemerisContext()
                # Randomly try different operations
                if thread_id % 3 == 0:
                    _ = ctx.get_loader()
                elif thread_id % 3 == 1:
                    _ = ctx.get_timescale()
                else:
                    _ = ctx.get_planets()
            except Exception as e:
                with errors_lock:
                    errors.append(f"Thread {thread_id}: {e}")

        def close_resources(thread_id: int):
            try:
                EphemerisContext.close()
            except Exception as e:
                with errors_lock:
                    errors.append(f"Close thread {thread_id}: {e}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads + 5) as ex:
            # Submit initialization tasks
            init_futures = [ex.submit(init_resources, i) for i in range(num_threads)]

            # Submit some close tasks in the middle
            close_futures = [ex.submit(close_resources, i) for i in range(5)]

            concurrent.futures.wait(init_futures + close_futures)

        # We expect no errors - the lock should prevent race conditions
        assert not errors, f"Errors during concurrent close/init: {errors}"


# =============================================================================
# TEST: get_planets() Ephemeris File Path Branches
# =============================================================================


class TestGetPlanetsEphemerisPathBranches:
    """
    Tests for get_planets() ephemeris file path search branches.

    The get_planets() method (context.py:146-180) has 4 code branches:
    1. Custom _SHARED_EPHE_PATH is set and file exists -> load from custom path
    2. Custom _SHARED_EPHE_PATH is set but file doesn't exist -> fall through
    3. File exists in data dir (~/.libephemeris/) -> load from data dir
    4. File doesn't exist anywhere -> download from internet (bare filename)

    These tests mock os.path.exists to control which branches are taken.
    """

    def test_branch_custom_ephe_path_file_exists(self, monkeypatch):
        """
        Branch 1: Custom _SHARED_EPHE_PATH is set and file exists.

        When a custom ephemeris path is configured and the file exists there,
        get_planets() should load from that custom path.
        """
        # Track which file path was loaded
        loaded_paths = []

        # Create a mock SpiceKernel
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        # Mock os.path.exists to return True for custom path
        custom_path = "/custom/ephemeris/path"
        custom_file_path = f"{custom_path}/de440.bsp"

        def mock_exists(path):
            return path == custom_file_path

        monkeypatch.setattr("os.path.exists", mock_exists)

        # Set up shared state with custom path
        ctx_module._SHARED_EPHE_PATH = custom_path
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        # Mock the loader to capture what path is loaded
        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        ctx = EphemerisContext()
        result = ctx.get_planets()

        # Verify the custom path was used
        assert len(loaded_paths) == 1
        assert loaded_paths[0] == custom_file_path
        assert result is mock_kernel

    def test_branch_custom_ephe_path_file_not_exists_fallback_to_data_dir(
        self, monkeypatch
    ):
        """
        Branch 2+3: Custom _SHARED_EPHE_PATH is set but file doesn't exist,
        then falls back to data dir (~/.libephemeris/) where file exists.

        Tests the fallback behavior when custom path is configured but
        the ephemeris file is not found there.
        """
        loaded_paths = []
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        custom_path = "/custom/ephemeris/path"
        custom_file_path = f"{custom_path}/de440.bsp"

        import os

        # Mock _get_data_dir to return a predictable path
        mock_data_dir = "/mock/libephemeris/data"
        data_dir_file_path = os.path.join(mock_data_dir, "de440.bsp")

        monkeypatch.setattr("libephemeris.state._get_data_dir", lambda: mock_data_dir)

        def mock_exists(path):
            # Custom path doesn't exist, but data dir does
            if path == custom_file_path:
                return False
            if path == data_dir_file_path:
                return True
            return False

        monkeypatch.setattr("os.path.exists", mock_exists)

        ctx_module._SHARED_EPHE_PATH = custom_path
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        ctx = EphemerisContext()
        result = ctx.get_planets()

        # Verify data dir path was used (fallback from custom path)
        assert len(loaded_paths) == 1
        assert loaded_paths[0] == data_dir_file_path
        assert result is mock_kernel

    def test_branch_no_custom_path_data_dir_exists(self, monkeypatch):
        """
        Branch 3: No custom path set, file exists in data dir.

        The default behavior when no custom ephemeris path is configured
        and the ephemeris file is found in the data dir (~/.libephemeris/).
        """
        loaded_paths = []
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        import os

        # Mock _get_data_dir to return a predictable path
        mock_data_dir = "/mock/libephemeris/data"
        data_dir_file_path = os.path.join(mock_data_dir, "de440.bsp")

        monkeypatch.setattr("libephemeris.state._get_data_dir", lambda: mock_data_dir)

        def mock_exists(path):
            return path == data_dir_file_path

        monkeypatch.setattr("os.path.exists", mock_exists)

        # No custom path set
        ctx_module._SHARED_EPHE_PATH = None
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        ctx = EphemerisContext()
        result = ctx.get_planets()

        # Verify data dir path was used
        assert len(loaded_paths) == 1
        assert loaded_paths[0] == data_dir_file_path
        assert result is mock_kernel

    def test_branch_no_file_anywhere_downloads_from_internet(self, monkeypatch):
        """
        Branch 4: File doesn't exist anywhere, download from internet.

        When the ephemeris file is not found in any local path, the loader
        should attempt to download it from the internet using just the filename.
        """
        loaded_paths = []
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        def mock_exists(path):
            # No file exists anywhere
            return False

        monkeypatch.setattr("os.path.exists", mock_exists)

        # No custom path set
        ctx_module._SHARED_EPHE_PATH = None
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        ctx = EphemerisContext()
        result = ctx.get_planets()

        # Verify just the filename was passed (triggering download)
        assert len(loaded_paths) == 1
        assert loaded_paths[0] == "de440.bsp"
        assert result is mock_kernel

    def test_custom_ephe_path_via_constructor(self, monkeypatch):
        """
        Test that ephe_path passed to constructor sets _SHARED_EPHE_PATH.

        Verifies that the constructor properly propagates the custom
        ephemeris path to the shared configuration.
        """
        loaded_paths = []
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        custom_path = "/user/provided/path"
        custom_file_path = f"{custom_path}/de440.bsp"

        def mock_exists(path):
            return path == custom_file_path

        monkeypatch.setattr("os.path.exists", mock_exists)

        # Start with no shared path
        ctx_module._SHARED_EPHE_PATH = None
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        # Create context with custom path
        ctx = EphemerisContext(ephe_path=custom_path)

        # Verify shared state was updated by constructor
        assert ctx_module._SHARED_EPHE_PATH == custom_path

        result = ctx.get_planets()

        # Verify custom path was used
        assert len(loaded_paths) == 1
        assert loaded_paths[0] == custom_file_path
        assert result is mock_kernel

    def test_custom_ephe_file_via_constructor(self, monkeypatch):
        """
        Test that ephe_file passed to constructor sets _SHARED_EPHE_FILE.

        Verifies that a custom ephemeris filename is properly used.
        The code checks data dir (~/.libephemeris/) for the file.
        """
        loaded_paths = []
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        import os

        custom_file = "de441.bsp"

        # Mock _get_data_dir to return a predictable path
        mock_data_dir = "/mock/libephemeris/data"
        data_dir_file_path = os.path.join(mock_data_dir, custom_file)

        monkeypatch.setattr("libephemeris.state._get_data_dir", lambda: mock_data_dir)

        def mock_exists(path):
            return path == data_dir_file_path

        monkeypatch.setattr("os.path.exists", mock_exists)

        ctx_module._SHARED_EPHE_PATH = None
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        # Create context with custom ephemeris file
        ctx = EphemerisContext(ephe_file=custom_file)

        # Verify shared state was updated by constructor
        assert ctx_module._SHARED_EPHE_FILE == custom_file

        result = ctx.get_planets()

        # Verify custom file was used via data dir
        assert len(loaded_paths) == 1
        assert loaded_paths[0] == data_dir_file_path
        assert result is mock_kernel

    def test_get_planets_returns_spicekernel(self):
        """
        Verify that get_planets() returns a SpiceKernel instance.

        This is an integration test that uses the real ephemeris file.
        """
        ctx = EphemerisContext()
        planets = ctx.get_planets()

        assert planets is not None
        assert isinstance(planets, SpiceKernel)

    def test_custom_ephe_path_with_custom_file(self, monkeypatch):
        """
        Test both custom path and custom file together.

        Verifies that both constructor parameters work in combination.
        """
        loaded_paths: list[str] = []
        mock_kernel = type("MockSpiceKernel", (), {})()

        def mock_loader_call(path):
            loaded_paths.append(path)
            return mock_kernel

        custom_path = "/opt/ephemeris"
        custom_file = "de441.bsp"
        custom_full_path = f"{custom_path}/{custom_file}"

        def mock_exists(path):
            return path == custom_full_path

        monkeypatch.setattr("os.path.exists", mock_exists)

        ctx_module._SHARED_EPHE_PATH = None
        ctx_module._SHARED_EPHE_FILE = "de440.bsp"

        mock_loader = type(
            "MockLoader", (), {"__call__": lambda s, p: mock_loader_call(p)}
        )()
        ctx_module._SHARED_LOADER = mock_loader

        # Create context with both custom path and file
        ctx = EphemerisContext(ephe_path=custom_path, ephe_file=custom_file)

        assert ctx_module._SHARED_EPHE_PATH == custom_path
        assert ctx_module._SHARED_EPHE_FILE == custom_file

        result = ctx.get_planets()

        assert len(loaded_paths) == 1
        assert loaded_paths[0] == custom_full_path
        assert result is mock_kernel


# =============================================================================
# TEST: SPK Body Registration (Context-local)
# =============================================================================


class TestContextSpkBodyRegistration:
    """
    Tests for EphemerisContext SPK body registration methods.

    These tests cover the register_spk_body(), unregister_spk_body(),
    get_spk_body_info(), and list_spk_bodies() methods defined in
    context.py:274-381.

    The tests use mocking to avoid network dependencies and ensure
    proper error handling for NAIF ID validation.
    """

    @pytest.fixture(autouse=True)
    def reset_spk_state(self):
        """Reset SPK state before and after each test."""
        from libephemeris import state

        # Save original state
        orig_body_map = dict(state._SPK_BODY_MAP)
        orig_kernels = dict(state._SPK_KERNELS)

        # Clear before test
        state._SPK_BODY_MAP.clear()
        state._SPK_KERNELS.clear()

        yield

        # Restore after test
        state._SPK_BODY_MAP.clear()
        state._SPK_BODY_MAP.update(orig_body_map)
        state._SPK_KERNELS.clear()
        state._SPK_KERNELS.update(orig_kernels)

    def test_register_spk_body_with_mock_file(self, tmp_path, monkeypatch):
        """
        Test registering a body with a mock SPK file.

        Verifies the complete registration flow:
        1. File path resolution
        2. Kernel loading (mocked)
        3. Context-local map update
        """
        from libephemeris import state

        # Create a mock SPK file
        mock_spk_path = tmp_path / "mock_body.bsp"
        mock_spk_path.touch()

        # Create mock kernel with a NAIF ID that exists
        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "mock_target"
                    if key == 2099999
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2099999"],
            },
        )()

        # Mock _load_spk_kernel to return our mock kernel
        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        # Register the body
        ctx = EphemerisContext()
        ctx.register_spk_body(
            ipl=15000,  # Custom body ID
            spk_file=str(mock_spk_path),
            naif_id=2099999,
        )

        # Verify registration in context-local map
        assert 15000 in ctx._spk_body_map
        assert ctx._spk_body_map[15000] == (str(mock_spk_path), 2099999)

    def test_get_spk_body_info_after_registration(self, tmp_path, monkeypatch):
        """
        Test get_spk_body_info() returns correct info after registration.

        Verifies that the method returns the correct (spk_file, naif_id) tuple.
        """
        from libephemeris import state

        mock_spk_path = tmp_path / "test_body.bsp"
        mock_spk_path.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2001234
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2001234"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()

        # Before registration
        assert ctx.get_spk_body_info(15001) is None

        # Register
        ctx.register_spk_body(
            ipl=15001,
            spk_file=str(mock_spk_path),
            naif_id=2001234,
        )

        # After registration
        info = ctx.get_spk_body_info(15001)
        assert info is not None
        assert info == (str(mock_spk_path), 2001234)

    def test_unregister_spk_body(self, tmp_path, monkeypatch):
        """
        Test unregister_spk_body() removes the registration.

        Verifies that after unregistration:
        1. The body is no longer in context-local map
        2. get_spk_body_info() returns None
        """
        from libephemeris import state

        mock_spk_path = tmp_path / "unregister_test.bsp"
        mock_spk_path.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2005678
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2005678"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()

        # Register
        ctx.register_spk_body(
            ipl=15002,
            spk_file=str(mock_spk_path),
            naif_id=2005678,
        )
        assert ctx.get_spk_body_info(15002) is not None

        # Unregister
        ctx.unregister_spk_body(15002)

        # Verify removal
        assert 15002 not in ctx._spk_body_map
        assert ctx.get_spk_body_info(15002) is None

    def test_unregister_nonexistent_body_no_error(self):
        """
        Test unregister_spk_body() on non-existent body does not raise.

        This tests the guard clause at context.py:344-345.
        """
        ctx = EphemerisContext()

        # Should not raise
        ctx.unregister_spk_body(99999)

    def test_list_spk_bodies_empty(self):
        """
        Test list_spk_bodies() returns empty dict when no bodies registered.
        """
        ctx = EphemerisContext()
        result = ctx.list_spk_bodies()
        assert result == {}

    def test_list_spk_bodies_with_context_local(self, tmp_path, monkeypatch):
        """
        Test list_spk_bodies() includes context-local registrations.
        """
        from libephemeris import state

        mock_spk_path = tmp_path / "list_test.bsp"
        mock_spk_path.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2009999
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2009999"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()
        ctx.register_spk_body(
            ipl=15003,
            spk_file=str(mock_spk_path),
            naif_id=2009999,
        )

        bodies = ctx.list_spk_bodies()
        assert 15003 in bodies
        assert bodies[15003] == (str(mock_spk_path), 2009999)

    def test_list_spk_bodies_merges_global_and_local(self, tmp_path, monkeypatch):
        """
        Test list_spk_bodies() merges global and context-local registrations.

        Context-local registrations should take precedence over global.
        """
        from libephemeris import state

        # Add a global registration
        state._SPK_BODY_MAP[10000] = ("/global/path.bsp", 2000001)

        mock_spk_path = tmp_path / "local.bsp"
        mock_spk_path.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2000002
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2000002"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()
        ctx.register_spk_body(
            ipl=10001,
            spk_file=str(mock_spk_path),
            naif_id=2000002,
        )

        bodies = ctx.list_spk_bodies()

        # Should include both global and local
        assert 10000 in bodies
        assert 10001 in bodies
        assert bodies[10000] == ("/global/path.bsp", 2000001)
        assert bodies[10001] == (str(mock_spk_path), 2000002)

    def test_context_local_overrides_global(self, tmp_path, monkeypatch):
        """
        Test that context-local registration overrides global for same ipl.
        """
        from libephemeris import state

        # Add global registration for ipl=10002
        state._SPK_BODY_MAP[10002] = ("/global/old.bsp", 2000100)

        mock_spk_path = tmp_path / "local_override.bsp"
        mock_spk_path.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2000200
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2000200"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()
        ctx.register_spk_body(
            ipl=10002,
            spk_file=str(mock_spk_path),
            naif_id=2000200,
        )

        # Context-local should override global
        info = ctx.get_spk_body_info(10002)
        assert info == (str(mock_spk_path), 2000200)

        # list_spk_bodies should show local, not global
        bodies = ctx.list_spk_bodies()
        assert bodies[10002] == (str(mock_spk_path), 2000200)

    def test_register_spk_body_file_not_found(self):
        """
        Test register_spk_body() raises FileNotFoundError for missing file.

        This tests the error path at context.py:308-312.
        """
        ctx = EphemerisContext()

        with pytest.raises(FileNotFoundError) as exc_info:
            ctx.register_spk_body(
                ipl=15004,
                spk_file="/nonexistent/path/to/file.bsp",
                naif_id=2000001,
            )

        assert "/nonexistent/path/to/file.bsp" in str(exc_info.value)

    def test_register_spk_body_invalid_naif_id(self, tmp_path, monkeypatch):
        """
        Test register_spk_body() raises ValueError for invalid NAIF ID.

        This tests the NAIF ID validation at context.py:319-332.
        When the NAIF ID is not found in the kernel, a ValueError should be
        raised with helpful information about available targets.
        """
        from libephemeris import state

        mock_spk_path = tmp_path / "invalid_naif.bsp"
        mock_spk_path.touch()

        # Create kernel that only has NAIF ID 2001000
        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (_ for _ in ()).throw(KeyError(key)),
                "names": lambda self: ["2001000", "2001001", "2001002"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()

        with pytest.raises(ValueError) as exc_info:
            ctx.register_spk_body(
                ipl=15005,
                spk_file=str(mock_spk_path),
                naif_id=9999999,  # Invalid NAIF ID
            )

        error_msg = str(exc_info.value)
        assert "9999999" in error_msg
        assert "not found" in error_msg.lower()

    def test_register_spk_body_naif_id_as_string(self, tmp_path, monkeypatch):
        """
        Test register_spk_body() accepts naif_id as string.

        This tests the string conversion at context.py:299-300.
        """
        from libephemeris import state

        mock_spk_path = tmp_path / "string_naif.bsp"
        mock_spk_path.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2003000
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2003000"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()
        ctx.register_spk_body(
            ipl=15006,
            spk_file=str(mock_spk_path),
            naif_id="2003000",  # String instead of int
        )

        # Should be stored as int
        info = ctx.get_spk_body_info(15006)
        assert info is not None
        assert info[1] == 2003000
        assert isinstance(info[1], int)

    def test_register_spk_body_relative_path_resolution(self, tmp_path, monkeypatch):
        """
        Test register_spk_body() resolves relative path via library path.

        This tests the path resolution logic at context.py:303-309.
        """
        from libephemeris import state

        # Create file in library path
        lib_path = tmp_path / "lib"
        lib_path.mkdir()
        spk_file = lib_path / "relative_test.bsp"
        spk_file.touch()

        # Mock get_library_path to return our temp dir
        monkeypatch.setattr(state, "get_library_path", lambda: str(lib_path))

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key == 2004000
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2004000"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx = EphemerisContext()
        ctx.register_spk_body(
            ipl=15007,
            spk_file="relative_test.bsp",  # Relative path
            naif_id=2004000,
        )

        # Should be registered with full path
        info = ctx.get_spk_body_info(15007)
        assert info is not None
        assert info[0] == str(spk_file)

    def test_context_isolation_between_instances(self, tmp_path, monkeypatch):
        """
        Test that SPK registrations are isolated between context instances.

        Each context should have its own _spk_body_map.
        """
        from libephemeris import state

        mock_spk_path1 = tmp_path / "ctx1.bsp"
        mock_spk_path1.touch()
        mock_spk_path2 = tmp_path / "ctx2.bsp"
        mock_spk_path2.touch()

        mock_kernel = type(
            "MockSpiceKernel",
            (),
            {
                "__getitem__": lambda self, key: (
                    "target"
                    if key in (2005000, 2006000)
                    else ((_ for _ in ()).throw(KeyError(key)))
                ),
                "names": lambda self: ["2005000", "2006000"],
            },
        )()

        def mock_load_spk_kernel(filepath):
            state._SPK_KERNELS[filepath] = mock_kernel
            return mock_kernel

        monkeypatch.setattr(state, "_load_spk_kernel", mock_load_spk_kernel)

        ctx1 = EphemerisContext()
        ctx2 = EphemerisContext()

        # Register different bodies in each context
        ctx1.register_spk_body(ipl=15008, spk_file=str(mock_spk_path1), naif_id=2005000)
        ctx2.register_spk_body(ipl=15009, spk_file=str(mock_spk_path2), naif_id=2006000)

        # ctx1 should only see 15008
        assert ctx1.get_spk_body_info(15008) is not None
        assert ctx1.get_spk_body_info(15009) is None

        # ctx2 should only see 15009
        assert ctx2.get_spk_body_info(15008) is None
        assert ctx2.get_spk_body_info(15009) is not None

    def test_get_spk_body_info_falls_back_to_global(self):
        """
        Test get_spk_body_info() falls back to global map.

        This tests the fallback at context.py:363-365.
        """
        from libephemeris import state

        # Add to global map only
        state._SPK_BODY_MAP[10003] = ("/global/fallback.bsp", 2000003)

        ctx = EphemerisContext()

        # Should find in global map
        info = ctx.get_spk_body_info(10003)
        assert info is not None
        assert info == ("/global/fallback.bsp", 2000003)
