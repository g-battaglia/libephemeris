# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Read original LEB2 bytes from bounded PostgreSQL blocks."""

from __future__ import annotations

import threading
from collections import OrderedDict

from libephemeris import CoefficientSourceError

from .config import RuntimeConfig
from .pool import get_pool

BLOCK_SIZE = 64 * 1024


class PgByteSource:
    """A read-at source for one immutable, content-addressed LEB2 file."""

    def __init__(
        self, name: str, sha256: str, size: int, config: RuntimeConfig
    ) -> None:
        if size <= 0:
            raise CoefficientSourceError("Invalid PostgreSQL file size")
        self.name = name
        self.locator = f"postgres://{sha256}"
        self._sha256 = sha256
        self._size = size
        self._config = config
        self._blocks: OrderedDict[int, bytes] = OrderedDict()
        self._lock = threading.Lock()
        self._closed = False

    def __len__(self) -> int:
        return self._size

    def read(self, offset: int, size: int) -> bytes:
        """Read an exact byte interval, batching cache misses in one query."""
        if offset < 0 or size < 0 or offset + size > self._size:
            raise CoefficientSourceError("Invalid PostgreSQL byte interval")
        if size == 0:
            return b""
        keys = range(offset // BLOCK_SIZE, (offset + size - 1) // BLOCK_SIZE + 1)
        with self._lock:
            if self._closed:
                raise CoefficientSourceError("PostgreSQL byte source is closed")
            blocks = {key: self._blocks[key] for key in keys if key in self._blocks}
            for key in blocks:
                self._blocks.move_to_end(key)
        missing = [key for key in keys if key not in blocks]
        if missing:
            try:
                with get_pool(self._config).connection() as conn:
                    rows = conn.execute(
                        "SELECT block_no, data FROM libephemeris.blocks "
                        "WHERE sha256=%s AND block_no=ANY(%s)",
                        (self._sha256, missing),
                    ).fetchall()
                fetched = {int(key): bytes(data) for key, data in rows}
            except Exception:
                raise CoefficientSourceError(
                    "PostgreSQL coefficient source unavailable"
                ) from None
            if set(fetched) != set(missing) or any(
                len(data) != min(BLOCK_SIZE, self._size - key * BLOCK_SIZE)
                for key, data in fetched.items()
            ):
                raise CoefficientSourceError("Missing or truncated PostgreSQL blocks")
            blocks.update(fetched)
            with self._lock:
                if not self._closed:
                    self._blocks.update(fetched)
                    while len(self._blocks) > self._config.cache_blocks:
                        self._blocks.popitem(last=False)
        data = b"".join(blocks[key] for key in keys)
        start = offset % BLOCK_SIZE
        return data[start : start + size]

    def close(self) -> None:
        """Release instance bytes without closing the shared process pool."""
        with self._lock:
            self._closed = True
            self._blocks.clear()
