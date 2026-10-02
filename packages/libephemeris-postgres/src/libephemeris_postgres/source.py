# SPDX-License-Identifier: AGPL-3.0-only
from __future__ import annotations

import threading

from libephemeris.leb2_reader import LEB2Reader

from .db import CoefficientSourceError, query

BLOCK_SIZE = 64 * 1024
CACHE_BLOCKS = 256


class PgByteSource:
    """Mmap-like slices of one published, content-addressed LEB2 file."""

    def __init__(self, sha256: str, size: int) -> None:
        if size <= 0:
            raise CoefficientSourceError("Invalid PostgreSQL file size")
        self._sha256, self._size = sha256, size
        self._blocks: dict[int, bytes] = {}
        self._lock = threading.Lock()
        self._closed = False

    def __len__(self) -> int:
        return self._size

    def __getitem__(self, interval: slice) -> bytes:
        """Batch missing blocks; retain returned bytes independently of cache eviction."""
        start, stop = interval.start, interval.stop
        if interval.step is not None or not 0 <= start <= stop <= self._size:
            raise CoefficientSourceError("Invalid PostgreSQL byte interval")
        with self._lock:
            if self._closed:
                raise CoefficientSourceError("PostgreSQL byte source is closed")
            if start == stop:
                return b""
            keys = range(start // BLOCK_SIZE, (stop - 1) // BLOCK_SIZE + 1)
            blocks = {key: self._blocks[key] for key in keys if key in self._blocks}
            missing = [key for key in keys if key not in blocks]
            if missing:
                fetched = dict(
                    query(
                        "SELECT block_no,data FROM libephemeris.blocks "
                        "WHERE sha256=%s AND block_no=ANY(%s)",
                        (self._sha256, missing),
                    )
                )
                if set(fetched) != set(missing) or any(
                    len(data) != min(BLOCK_SIZE, self._size - key * BLOCK_SIZE)
                    for key, data in fetched.items()
                ):
                    raise CoefficientSourceError(
                        "Missing or truncated PostgreSQL blocks"
                    )
                blocks.update(fetched)
                if len(self._blocks) + len(fetched) > CACHE_BLOCKS:
                    self._blocks.clear()
                self._blocks.update(list(fetched.items())[-CACHE_BLOCKS:])
            return b"".join(blocks[key] for key in keys)[
                start % BLOCK_SIZE : start % BLOCK_SIZE + stop - start
            ]

    def close(self) -> None:
        """Release instance bytes, not the process connection."""
        with self._lock:
            self._closed = True
            self._blocks.clear()


class PgLEB2Reader(LEB2Reader):
    """Keep damaged remote chunks fatal instead of entering scientific fallback."""

    def _decompress_chunk_uncached(self, body_id: int, chunk_idx: int) -> bytes:
        try:
            return super()._decompress_chunk_uncached(body_id, chunk_idx)
        except Exception:
            raise CoefficientSourceError(
                "Could not decode PostgreSQL LEB2 source"
            ) from None
