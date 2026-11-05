"""Bounded asyncio offload to the native engine; Python 3.12 or later.

Run with an installed v4 wheel:
    python python/examples/async_worker.py /path/to/leb/base_core.leb2 --workers 4

Initialize resources before serving. Cancelling an awaiting request leaves its
submitted calculation and admission slot owned until the worker really ends.
Shutdown must finish before the application's event loop closes.
"""
from __future__ import annotations

import argparse
import asyncio
from pathlib import Path

import libephemeris as ephe
from bounded_workers import BoundedWorkers


class AsyncEphemerisWorkers:
    """Reuse one thread-safe native context with bounded submitted work."""

    def __init__(self, context: ephe.EphemerisContext, workers: int = 1):
        self._context = context
        self._workers = BoundedWorkers(workers=workers, max_queued=2 * workers)

    async def calculate(self, jd_ut1: float, body: int, flags: int = 258):
        """Await a UT1 calculation without running it on the event loop."""
        return await self._workers.run(self._context.calc_ut, jd_ut1, body, flags)

    async def close(self):
        """Stop admission and drain submitted calculations without blocking asyncio."""
        await self._workers.close()

    async def __aenter__(self):
        return self

    async def __aexit__(self, *_):
        await self.close()


def initialized_context(leb_file: Path):
    """Configure once before starting the event loop or accepting requests."""
    ephe.set_calc_mode("leb")
    ephe.set_precision_tier("base")
    ephe.set_network_policy("sealed")
    ephe.set_iers_delta_t_enabled(False)
    context = ephe.EphemerisContext()
    context.set_leb_file(str(leb_file.resolve()))
    context.calc_ut(2451545.0, ephe.SUN, ephe.FLG_SWIEPH | ephe.FLG_SPEED)
    return context


async def main(context, workers):
    async with AsyncEphemerisWorkers(context, workers) as pool:
        # The example stays below admission capacity. An HTTP application maps
        # ServerBusyError to 503 instead of accumulating waiting requests.
        results = await asyncio.gather(*(
            pool.calculate(2451545.0 + day, ephe.MOON) for day in range(workers)
        ))
        print(results[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("leb_file", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    arguments = parser.parse_args()
    if arguments.workers < 1:
        parser.error("--workers must be positive")
    asyncio.run(main(initialized_context(arguments.leb_file), arguments.workers))
