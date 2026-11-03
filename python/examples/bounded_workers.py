"""Bound whole synchronous jobs outside an asyncio application's event loop.

This application adapter is not part of the ephemeris runtime package. One
admitted job owns its capacity until the underlying thread really completes.
Full admission fails immediately, so overload cannot grow a queue of waiters.
"""
from __future__ import annotations

import asyncio
import contextvars
import functools
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import ParamSpec, TypeVar

P = ParamSpec("P")
T = TypeVar("T")


class ServerBusyError(RuntimeError):
    """The configured running-plus-queued calculation capacity is full."""


class BoundedWorkers:
    """A process-local pool with immediate, strictly bounded admission.

    Create, use and close on one event loop. Submit a complete subject/chart
    operation, keeping configuration, native trace scopes and cleanup on the
    same worker thread. ContextVars are copied; native thread-local tokens are
    deliberately not transported. Call ``close`` before the event loop ends.
    """

    def __init__(self, workers: int = 1, max_queued: int = 2):
        if workers < 1 or max_queued < 0:
            raise ValueError("workers must be positive and max_queued nonnegative")
        self._loop = asyncio.get_running_loop()
        self._pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="ephemeris")
        self._capacity = workers + max_queued
        self._pending = 0
        self._closed = False
        self._idle = asyncio.Event()
        self._idle.set()

    async def run(self, function: Callable[P, T], /, *args: P.args, **kwargs: P.kwargs) -> T:
        """Run a whole job; reject overload without waiting for admission."""
        if asyncio.get_running_loop() is not self._loop:
            raise RuntimeError("worker pool belongs to a different event loop")
        if self._closed:
            raise RuntimeError("worker pool is closed")
        if self._pending >= self._capacity:
            raise ServerBusyError("calculation capacity is full; retry later")
        # No await separates the check and reservation on the owning loop.
        self._pending += 1
        self._idle.clear()
        try:
            call = functools.partial(contextvars.copy_context().run, function, *args, **kwargs)
            result = self._loop.create_future()
            work = self._pool.submit(call)
        except BaseException:
            self._finished()
            raise
        # Keep completion ownership separate from the request's asyncio future.
        # Cancelling that future cannot cancel a running thread or free capacity.
        work.add_done_callback(
            lambda future: self._loop.call_soon_threadsafe(self._complete, future, result)
        )
        try:
            return await result
        except asyncio.CancelledError:
            # This cancels a queued job if it has not started. Running native
            # work continues and retains its slot until actual completion.
            work.cancel()
            raise

    def _finished(self) -> None:
        self._pending -= 1
        if self._pending == 0:
            self._idle.set()

    def _complete(self, work, result: asyncio.Future) -> None:
        try:
            if work.cancelled():
                result.cancel()
                return
            # Retrieve late failures even after a client disconnects. The live
            # caller receives the original exception; a cancelled one has no
            # orphaned asyncio future carrying request data into default logs.
            error = work.exception()
            if not result.done():
                if error is None:
                    result.set_result(work.result())
                else:
                    result.set_exception(error)
        finally:
            self._finished()

    async def close(self) -> None:
        """Stop admission and drain actual jobs before closing the executor."""
        self._closed = True
        await self._idle.wait()
        await asyncio.to_thread(self._pool.shutdown, wait=True, cancel_futures=False)

    async def __aenter__(self):
        return self

    async def __aexit__(self, *_):
        await self.close()
