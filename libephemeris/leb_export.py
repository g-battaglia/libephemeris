# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Read-only export of decoded native LEB coefficient segments.

Provenance:
    Project-native binary layouts and chunk indexes are defined in
    ``leb_format``. Export preserves the decoded IEEE-754 bytes, without
    fitting or changing coefficients. Chunk selection uses integer indexes.
"""

from __future__ import annotations

import struct
from collections.abc import Iterator
from dataclasses import replace
from typing import Any

from .exceptions import LEBCorruptionError
from .leb_format import BodyEntry, CompressedBodyEntry, NutationHeader, StarEntry


class LEBExportMixin:
    """Public coefficient export methods shared by native file readers.

    These methods are independent of the readers' evaluation hot paths.
    Body id -1 denotes the auxiliary nutation series in page exports only.
    """

    # Native readers own these fields. Export never participates in evaluation.
    _mm: Any
    _bodies: Any
    _nutation: NutationHeader | None
    _nutation_data_offset: int
    _delta_t_jds: list[float]
    _delta_t_vals: list[float]
    _stars: dict[int, StarEntry]
    _chunked: bool
    _chunk_index: Any
    _decompress_chunk: Any
    _decompress_body: Any
    jd_range: Any

    @property
    def header_jd_range(self) -> tuple[float, float]:
        """Return the original artifact header interval."""
        return self.jd_range

    @property
    def bodies(self) -> dict[int, BodyEntry | CompressedBodyEntry]:
        """Return detached copies of the native body metadata."""
        return {key: replace(value) for key, value in self._bodies.items()}

    @property
    def nutation_header(self) -> NutationHeader | None:
        """Return detached nutation metadata, if present."""
        return replace(self._nutation) if self._nutation is not None else None

    @property
    def delta_t_table(self) -> tuple[tuple[float, ...], tuple[float, ...]]:
        """Return immutable Julian Day and Delta-T (days) tables."""
        return tuple(self._delta_t_jds), tuple(self._delta_t_vals)

    @property
    def stars(self) -> dict[int, StarEntry]:
        """Return detached copies of the native star catalog."""
        return {key: replace(value) for key, value in self._stars.items()}

    def _export_segment_bytes(self, body_id: int, idx: int) -> bytes:
        mm = self._mm
        if mm is None:
            raise ValueError("LEB reader is closed")
        if body_id == -1:
            entry = self._nutation
            if entry is None:
                raise KeyError("No nutation data in this LEB file")
            offset = self._nutation_data_offset
        else:
            entry = self._bodies[body_id]
            offset = entry.data_offset
        if not isinstance(idx, int) or not 0 <= idx < entry.segment_count:
            raise IndexError("Coefficient segment index outside stored range")
        size = (entry.degree + 1) * entry.components * 8
        if body_id != -1 and isinstance(entry, CompressedBodyEntry):
            if self._chunked:
                chunks = self._chunk_index[body_id]
                lo, hi = 0, len(chunks)
                while lo < hi:
                    mid = (lo + hi) // 2
                    if chunks[mid].segment_start <= idx:
                        lo = mid + 1
                    else:
                        hi = mid
                chunk_idx = lo - 1
                if chunk_idx < 0:
                    raise LEBCorruptionError("Missing LEB2 segment chunk")
                chunk = chunks[chunk_idx]
                local = idx - chunk.segment_start
                if local >= chunk.segment_count:
                    raise LEBCorruptionError("Missing LEB2 segment chunk")
                data = self._decompress_chunk(body_id, chunk_idx)
                offset = local * size
            else:
                data = self._decompress_body(body_id)
                offset = idx * size
        else:
            data = mm
            offset += idx * size
        if offset < 0 or offset + size > len(data):
            raise LEBCorruptionError("Truncated exported coefficient segment")
        return bytes(data[offset : offset + size])

    def segment_coefficients(self, body_id: int, idx: int) -> tuple[float, ...]:
        """Export a body's segment by global integer index."""
        if body_id == -1:
            raise KeyError("Use nutation_coefficients for nutation")
        data = self._export_segment_bytes(body_id, idx)
        return struct.unpack(f"<{len(data) // 8}d", data)

    def nutation_coefficients(self, idx: int) -> tuple[float, ...]:
        """Export one nutation segment in component-major order."""
        data = self._export_segment_bytes(-1, idx)
        return struct.unpack(f"<{len(data) // 8}d", data)

    def iter_segment_pages(
        self, body_id: int, page_size: int = 32
    ) -> Iterator[tuple[int, int, bytes]]:
        """Yield (page number, first segment index, little-endian bytes).

        Args:
            body_id: Native body id, or -1 for nutation.
            page_size: Maximum number of consecutive segments per page.
        """
        if not isinstance(page_size, int) or page_size <= 0:
            raise ValueError("page_size must be a positive integer")
        entry = self._nutation if body_id == -1 else self._bodies[body_id]
        if entry is None:
            raise KeyError("No nutation data in this LEB file")
        for first in range(0, entry.segment_count, page_size):
            last = min(first + page_size, entry.segment_count)
            yield (
                first // page_size,
                first,
                b"".join(
                    self._export_segment_bytes(body_id, idx)
                    for idx in range(first, last)
                ),
            )
