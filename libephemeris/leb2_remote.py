# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Reuse native LEB2 evaluation over original remote file bytes."""

from __future__ import annotations

from typing import Protocol

from .exceptions import CoefficientSourceError
from .leb2_reader import LEB2Reader


class ByteSource(Protocol):
    """Read-at contract; close releases only instance resources."""

    name: str
    locator: str

    def __len__(self) -> int: ...
    def read(self, offset: int, size: int) -> bytes: ...
    def close(self) -> None: ...


class RemoteLEB2Reader(LEB2Reader):
    """Parse and evaluate provider-owned bytes with the native reader."""

    def __init__(self, source: ByteSource, *, reviewed: bool) -> None:
        self._source = source
        try:
            self._initialize(
                f"{source.locator.rstrip('/')}/{source.name}", source, None
            )
            if not self._chunked:
                raise ValueError("Remote sources require LEB2 v2")
            self._manifest_verified = reviewed
        except Exception:
            self.close()
            raise CoefficientSourceError("Could not open remote LEB2 source") from None

    def _read_blob(self, offset: int, size: int, what: str) -> bytes:
        try:
            if (
                self._mm is None
                or offset < 0
                or size < 0
                or offset + size > len(self._source)
            ):
                raise ValueError("Invalid byte interval")
            data = self._source.read(offset, size)
            if not isinstance(data, bytes) or len(data) != size:
                raise ValueError("Short byte read")
            return data
        except Exception:
            raise CoefficientSourceError("Could not read remote LEB2 source") from None

    def _decompress_chunk_uncached(self, body_id: int, chunk_idx: int) -> bytes:
        try:
            return super()._decompress_chunk_uncached(body_id, chunk_idx)
        except Exception:
            raise CoefficientSourceError(
                "Could not decode remote LEB2 source"
            ) from None

    def warm(self, jd_start: float, jd_end: float) -> None:
        """Remote bytes are fetched lazily, not mmap-prefaulted."""

    def cool(self) -> None:
        """The provider owns its bounded byte cache."""
