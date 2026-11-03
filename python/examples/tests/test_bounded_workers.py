"""Admission, cancellation and isolation of complete asynchronous worker jobs."""
from __future__ import annotations

import asyncio
import contextvars
import gc
import importlib.util
from pathlib import Path
import threading
import weakref

import pytest

_source = Path(__file__).resolve().parents[3] / "python/examples/bounded_workers.py"
_spec = importlib.util.spec_from_file_location("bounded_workers", _source)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
BoundedWorkers = _module.BoundedWorkers
ServerBusyError = _module.ServerBusyError


def test_overload_is_immediate_and_running_cancellation_keeps_capacity():
    async def scenario():
        started = asyncio.Event()
        release = threading.Event()
        loop = asyncio.get_running_loop()

        def blocked():
            loop.call_soon_threadsafe(started.set)
            assert release.wait(5)

        pool = BoundedWorkers(workers=1, max_queued=0)
        try:
            request = asyncio.create_task(pool.run(blocked))
            await asyncio.wait_for(started.wait(), 2)
            with pytest.raises(ServerBusyError):
                await pool.run(lambda: None)
            request.cancel()
            with pytest.raises(asyncio.CancelledError):
                await request
            # A disconnected client cannot make room for a second native job.
            with pytest.raises(ServerBusyError):
                await pool.run(lambda: None)
            # The loop still executes while the synchronous worker is blocked.
            tick = asyncio.Event()
            loop.call_soon(tick.set)
            await asyncio.wait_for(tick.wait(), 1)
        finally:
            release.set()
            await pool.close()
        assert pool._pending == 0

    asyncio.run(scenario())


def test_cancelled_queued_job_never_runs_and_shutdown_drains_running_job():
    async def scenario():
        started = asyncio.Event()
        release = threading.Event()
        queued_ran = threading.Event()
        loop = asyncio.get_running_loop()

        def blocked():
            loop.call_soon_threadsafe(started.set)
            assert release.wait(5)
            return "complete"

        pool = BoundedWorkers(workers=1, max_queued=1)
        try:
            running = asyncio.create_task(pool.run(blocked))
            await asyncio.wait_for(started.wait(), 2)
            queued = asyncio.create_task(pool.run(queued_ran.set))
            await asyncio.sleep(0)
            with pytest.raises(ServerBusyError):
                await pool.run(lambda: None)
            queued.cancel()
            with pytest.raises(asyncio.CancelledError):
                await queued
            shutdown = asyncio.create_task(pool.close())
            await asyncio.sleep(0)
            assert not shutdown.done()
            with pytest.raises(RuntimeError, match="closed"):
                await pool.run(lambda: None)
            release.set()
            assert await running == "complete"
            await shutdown
            assert not queued_ran.is_set()
        finally:
            release.set()
            await pool.close()

    asyncio.run(scenario())


def test_contexts_and_errors_are_isolated_and_recover_capacity():
    marker = contextvars.ContextVar("request_marker", default="outside")

    async def scenario():
        async with BoundedWorkers(workers=2, max_queued=0) as pool:
            async def request(value):
                token = marker.set(value)
                try:
                    return await pool.run(marker.get)
                finally:
                    marker.reset(token)

            assert await asyncio.gather(request("first"), request("second")) == ["first", "second"]
            assert marker.get() == "outside"

            def fail():
                raise ValueError("worker failure")

            with pytest.raises(ValueError, match="worker failure"):
                await pool.run(fail)
            assert await pool.run(lambda: "recovered") == "recovered"

    asyncio.run(scenario())


def test_failure_after_client_cancellation_is_observed():
    async def scenario():
        started = asyncio.Event()
        release = threading.Event()
        loop = asyncio.get_running_loop()
        unhandled = []
        previous = loop.get_exception_handler()
        loop.set_exception_handler(lambda _, context: unhandled.append(context))

        def fail_later():
            loop.call_soon_threadsafe(started.set)
            assert release.wait(5)
            raise ValueError("late worker failure")

        try:
            async with BoundedWorkers(workers=1, max_queued=0) as pool:
                request = asyncio.create_task(pool.run(fail_later))
                await asyncio.wait_for(started.wait(), 2)
                request.cancel()
                with pytest.raises(asyncio.CancelledError):
                    await request
                release.set()
            await asyncio.sleep(0)
            assert not unhandled
        finally:
            release.set()
            loop.set_exception_handler(previous)

    asyncio.run(scenario())


def test_completed_jobs_release_request_contexts():
    class RequestState:
        pass

    marker = contextvars.ContextVar("request_state")

    async def scenario():
        references = []
        async with BoundedWorkers(workers=2, max_queued=2) as pool:
            for _ in range(128):
                value = RequestState()
                references.append(weakref.ref(value))
                token = marker.set(value)
                try:
                    assert await pool.run(lambda: isinstance(marker.get(), RequestState))
                finally:
                    marker.reset(token)
                del token, value
            await asyncio.sleep(0)
            gc.collect()
            assert all(reference() is None for reference in references)

    asyncio.run(scenario())
