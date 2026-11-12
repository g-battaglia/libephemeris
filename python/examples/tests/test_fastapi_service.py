"""FastAPI readiness, worker reuse and startup cleanup without astronomical data."""
from __future__ import annotations

import asyncio
from contextlib import suppress
import importlib.util
from pathlib import Path
import sys
import threading

import pytest


ROOT = Path(__file__).resolve().parents[3]


@pytest.fixture
def app_and_pools(monkeypatch):
    def load(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        monkeypatch.setitem(sys.modules, name, module)
        spec.loader.exec_module(module)
        return module

    workers = load("bounded_workers", ROOT / "python/examples/bounded_workers.py")
    example = load("_fastapi_lifespan_example", ROOT / "python/examples/fastapi_service.py")
    pools = []

    class ObservedWorkers(workers.BoundedWorkers):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            pools.append(self)

    monkeypatch.setattr(example, "BoundedWorkers", ObservedWorkers)
    monkeypatch.setenv("EPHEMERIS_WORKERS", "1")
    monkeypatch.setenv("EPHEMERIS_MAX_QUEUED", "0")
    return example, pools


def test_readiness_waits_and_startup_reuses_the_only_calculation_worker(app_and_pools, monkeypatch):
    example, pools = app_and_pools

    async def scenario():
        loop = asyncio.get_running_loop()
        started = asyncio.Event()
        ready = asyncio.Event()
        stop = asyncio.Event()
        release = threading.Event()
        worker_ids = []
        offloads = []
        original_to_thread = asyncio.to_thread

        async def observe_to_thread(function, *args, **kwargs):
            offloads.append(function)
            return await original_to_thread(function, *args, **kwargs)

        def initialize():
            worker_ids.append(threading.get_ident())
            assert threading.current_thread().name.startswith("ephemeris")
            loop.call_soon_threadsafe(started.set)
            assert release.wait(5)

        async def serve():
            async with example.app.router.lifespan_context(example.app):
                worker_ids.append(await example.app.state.calculations.run(threading.get_ident))
                ready.set()
                await stop.wait()

        monkeypatch.setattr(example, "initialize", initialize)
        monkeypatch.setattr(asyncio, "to_thread", observe_to_thread)
        task = asyncio.create_task(serve())
        try:
            await asyncio.wait_for(started.wait(), 2)
            tick = asyncio.Event()
            loop.call_soon(tick.set)
            await asyncio.wait_for(tick.wait(), 1)
            assert not ready.is_set()
            assert not hasattr(example.app.state, "calculations")
            assert not offloads
            release.set()
            await asyncio.wait_for(ready.wait(), 2)
            assert len(worker_ids) == 2
            assert worker_ids[0] == worker_ids[1] != threading.get_ident()
            assert not offloads
        finally:
            release.set()
            stop.set()
            await asyncio.wait_for(task, 5)
        pool, = pools
        assert offloads == [pool._pool.shutdown]
        assert pool._closed and pool._pending == 0
        assert all(not thread.is_alive() for thread in pool._pool._threads)

    asyncio.run(scenario())


@pytest.mark.parametrize("cancel_startup", [False, True], ids=["failure", "cancellation"])
def test_failed_or_cancelled_startup_drains_and_closes(app_and_pools, monkeypatch, cancel_startup):
    example, pools = app_and_pools

    class StartupFailure(RuntimeError):
        pass

    async def scenario():
        loop = asyncio.get_running_loop()
        started = asyncio.Event()
        release = threading.Event()
        ready = False
        error = StartupFailure("initialization failed")
        unhandled = []
        previous_handler = loop.get_exception_handler()
        loop.set_exception_handler(lambda _, context: unhandled.append(context))

        def initialize():
            loop.call_soon_threadsafe(started.set)
            assert release.wait(5)
            raise error

        async def serve():
            nonlocal ready
            async with example.app.router.lifespan_context(example.app):
                ready = True

        monkeypatch.setattr(example, "initialize", initialize)
        task = asyncio.create_task(serve())
        try:
            await asyncio.wait_for(started.wait(), 2)
            pool, = pools
            assert pool._pending == 1
            if cancel_startup:
                task.cancel()
                await asyncio.sleep(0)
                assert pool._closed
                assert not task.done()
                assert pool._pending == 1
                release.set()
                with pytest.raises(asyncio.CancelledError):
                    await asyncio.wait_for(task, 2)
            else:
                release.set()
                with pytest.raises(StartupFailure) as caught:
                    await asyncio.wait_for(task, 2)
                assert caught.value is error
            await asyncio.sleep(0)
            assert not ready
            assert not hasattr(example.app.state, "calculations")
            assert pool._closed and pool._pending == 0
            assert all(not thread.is_alive() for thread in pool._pool._threads)
            assert not unhandled
        finally:
            release.set()
            with suppress(StartupFailure, asyncio.CancelledError):
                await asyncio.wait_for(task, 5)
            loop.set_exception_handler(previous_handler)

    asyncio.run(scenario())


@pytest.mark.parametrize("workers,queued", [("0", "0"), ("1", "-1"),
                                           ("invalid", "0"), ("1", "invalid")])
def test_invalid_configuration_is_rejected_before_data_acquisition(app_and_pools, monkeypatch,
                                                                  workers, queued):
    example, pools = app_and_pools
    acquisitions = []
    monkeypatch.setenv("EPHEMERIS_WORKERS", workers)
    monkeypatch.setenv("EPHEMERIS_MAX_QUEUED", queued)
    monkeypatch.setattr(example, "initialize", lambda: acquisitions.append("initialize"))

    async def scenario():
        with pytest.raises(ValueError):
            async with example.app.router.lifespan_context(example.app):
                pytest.fail("invalid startup reached readiness")
        assert not acquisitions
        assert not pools
        assert not hasattr(example.app.state, "calculations")

    asyncio.run(scenario())
