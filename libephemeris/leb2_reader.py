# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Reader for LEB2 compressed ephemeris files (v1 and v2).

Provides the same interface as LEBReader.

- **v1 (legacy)**: Decompresses each body's entire coefficient blob lazily on
  first access and caches the result.  Fast for small bodies but can be slow
  for large ones (e.g. Moon at 307 MB decompressed).

- **v2 (chunked)**: Each body is split into 10-year temporal chunks, each
  compressed independently.  Only the chunk containing the requested JD is
  decompressed (~300 KB for Moon), giving near-instant cold-start performance.
"""

from __future__ import annotations

import mmap
import os
import struct
import threading
from bisect import bisect_right
from typing import Dict, List, Optional, Tuple

from .leb_compression import decompress_body
from .leb_format import (
    CHUNK_ENTRY_SIZE,
    CHUNK_INDEX_HEADER_SIZE,
    COMPRESSED_BODY_ENTRY_SIZE,
    COORD_ECLIPTIC,
    COORD_GEO_ECLIPTIC,
    COORD_HELIO_ECL,
    DELTA_T_ENTRY_FMT,
    DELTA_T_ENTRY_SIZE,
    DELTA_T_HEADER_FMT,
    DELTA_T_HEADER_SIZE,
    HEADER_SIZE,
    LEB2_MAGIC,
    LEB2_VERSION,
    LEB2_VERSION_V1,
    NUTATION_HEADER_SIZE,
    SECTION_BODY_INDEX,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    SECTION_STARS,
    STAR_ENTRY_SIZE,
    ChunkEntry,
    CompressedBodyEntry,
    NutationHeader,
    SectionEntry,
    StarEntry,
    read_chunk_entry,
    read_chunk_index_header,
    read_compressed_body_entry,
    read_header,
    read_nutation_header,
    read_section_dir,
    read_star_entry,
    LEBCorruptionError,
    validate_entry_count,
    _madvise_ranges,
    _madvise_dontneed,
)
from .leb_reader import _clenshaw, _clenshaw_with_derivative


class LEB2Reader:
    """Reader for LEB2 compressed .leb2 files (v1 monolithic and v2 chunked).

    Same interface as LEBReader.

    Usage:
        with LEB2Reader("ephemeris.leb2") as reader:
            pos, vel = reader.eval_body(SUN, jd_tt)
    """

    def __init__(self, path: str) -> None:
        if not os.path.exists(path):
            raise FileNotFoundError(f"LEB file not found: {path}")

        self._path = path
        self._file = open(path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._cache: Dict[int, bytes] = {}  # v1: body_id -> full decompressed data
        self._chunk_cache: Dict[Tuple[int, int], bytes] = {}  # v2: (body_id, chunk_idx) -> decompressed chunk
        # Guards _cache/_chunk_cache and serializes decompression: chunks
        # decompress in ~0.05 ms, so a plain double-checked lock is adequate
        # and concurrent requests for the same chunk never redo the zstd
        # work. (_eval_cache stays lock-free: a tiny memo of result tuples
        # whose dict get/set are atomic under the GIL.)
        self._decomp_lock = threading.Lock()
        self._chunk_index: Dict[int, List[ChunkEntry]] = {}  # v2: body_id -> list of ChunkEntry
        self._chunked: bool = False  # True for v2 chunked format
        self._eval_cache: Dict[
            Tuple[int, float],
            Tuple[Tuple[float, float, float], Tuple[float, float, float]],
        ] = {}

        try:
            self._parse()
        except (struct.error, IndexError, KeyError) as exc:
            # Normalize truncated/corrupted-file failures to ValueError so
            # validation and the Skyfield fallback see one failure type.
            self.close()
            raise LEBCorruptionError(f"Corrupted LEB2 file {path!r}: {exc}") from exc
        except (OSError, ValueError):
            self.close()
            raise

    def _parse(self) -> None:
        self._header = read_header(self._mm, 0)
        if self._header.magic != LEB2_MAGIC:
            raise ValueError(
                f"Invalid LEB2 magic: {self._header.magic!r} (expected {LEB2_MAGIC!r})"
            )
        if self._header.version not in (LEB2_VERSION_V1, LEB2_VERSION):
            raise ValueError(
                f"Unsupported LEB2 version: {self._header.version} "
                f"(supported: {LEB2_VERSION_V1}, {LEB2_VERSION})"
            )
        self._chunked = self._header.version >= LEB2_VERSION

        # The header flags carry the compression algorithm ID; feeding a
        # COMPRESSION_NONE (or unknown-scheme) file to the zstd pipeline
        # would fail confusingly deep in decompression.
        from .leb_format import COMPRESSION_ZSTD_TRUNC_SHUFFLE

        if self._header.flags != COMPRESSION_ZSTD_TRUNC_SHUFFLE:
            raise ValueError(
                f"Unsupported LEB2 compression scheme: {self._header.flags}"
            )

        # Parse section directory
        self._sections: Dict[int, SectionEntry] = {}
        for i in range(self._header.section_count):
            offset = HEADER_SIZE + i * SECTION_DIR_SIZE
            sec = read_section_dir(self._mm, offset)
            self._sections[sec.section_id] = sec

        # Parse body index (CompressedBodyEntry, 68 bytes each)
        self._bodies: Dict[int, CompressedBodyEntry] = {}
        if SECTION_BODY_INDEX in self._sections:
            sec = self._sections[SECTION_BODY_INDEX]
            validate_entry_count(
                self._header.body_count,
                COMPRESSED_BODY_ENTRY_SIZE,
                sec.size,
                "LEB2 body index",
            )
            for i in range(self._header.body_count):
                offset = sec.offset + i * COMPRESSED_BODY_ENTRY_SIZE
                entry = read_compressed_body_entry(self._mm, offset)
                self._bodies[entry.body_id] = entry

        # For v2: parse chunk indices for each body
        if self._chunked:
            for body_id, entry in self._bodies.items():
                self._chunk_index[body_id] = self._parse_chunk_index(entry)

        # Parse nutation header (uncompressed, same as LEB1)
        self._nutation: Optional[NutationHeader] = None
        self._nutation_data_offset: int = 0
        if SECTION_NUTATION in self._sections:
            sec = self._sections[SECTION_NUTATION]
            self._nutation = read_nutation_header(self._mm, sec.offset)
            self._nutation_data_offset = sec.offset + NUTATION_HEADER_SIZE

        # Parse Delta-T table (uncompressed, same as LEB1)
        self._delta_t_jds: List[float] = []
        self._delta_t_vals: List[float] = []
        if SECTION_DELTA_T in self._sections:
            sec = self._sections[SECTION_DELTA_T]
            n_entries, _ = struct.unpack_from(DELTA_T_HEADER_FMT, self._mm, sec.offset)
            data_offset = sec.offset + DELTA_T_HEADER_SIZE
            for i in range(n_entries):
                off = data_offset + i * DELTA_T_ENTRY_SIZE
                jd, dt = struct.unpack_from(DELTA_T_ENTRY_FMT, self._mm, off)
                self._delta_t_jds.append(jd)
                self._delta_t_vals.append(dt)

        # Parse star catalog (uncompressed, same as LEB1)
        self._stars: Dict[int, StarEntry] = {}
        if SECTION_STARS in self._sections:
            sec = self._sections[SECTION_STARS]
            n_stars = sec.size // STAR_ENTRY_SIZE
            for i in range(n_stars):
                offset = sec.offset + i * STAR_ENTRY_SIZE
                star = read_star_entry(self._mm, offset)
                self._stars[star.star_id] = star

    def _parse_chunk_index(self, entry: CompressedBodyEntry) -> List[ChunkEntry]:
        """Parse the chunk index for a v2 body."""
        offset = entry.data_offset
        chunk_count, _ = read_chunk_index_header(self._mm, offset)
        offset += CHUNK_INDEX_HEADER_SIZE

        # Corrupted-entry guard: an inflated u32 chunk_count (up to ~4e9)
        # would otherwise allocate billions of ChunkEntry objects and read
        # far past the chunk index.
        validate_entry_count(
            chunk_count,
            CHUNK_ENTRY_SIZE,
            len(self._mm) - offset,
            f"LEB2 chunk index (body {entry.body_id})",
        )

        chunks = []
        for i in range(chunk_count):
            chunk = read_chunk_entry(self._mm, offset + i * CHUNK_ENTRY_SIZE)
            chunks.append(chunk)
        return chunks

    def _find_chunk(self, body_id: int, jd: float) -> int:
        """Find the chunk index containing jd via binary search."""
        chunks = self._chunk_index[body_id]
        # Binary search: chunks are sorted by jd_start
        lo, hi = 0, len(chunks) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if jd >= chunks[mid].jd_end:
                lo = mid + 1
            else:
                hi = mid
        return lo

    def _read_blob(self, offset: int, size: int, what: str) -> bytes:
        """Slice `size` bytes from the mmap, rejecting truncated files.

        An mmap slice past EOF silently returns fewer bytes — surface a
        clear truncation error instead of an opaque zstd failure.

        Args:
            offset: Byte offset of the blob in the file.
            size: Expected blob size in bytes.
            what: Human-readable label for the error message.

        Raises:
            LEBCorruptionError: If fewer than `size` bytes are available.
        """
        blob = self._mm[offset : offset + size]
        if len(blob) < size:
            raise LEBCorruptionError(
                f"Truncated LEB2 file: {what} needs {size} bytes at "
                f"offset {offset}, only {len(blob)} available"
            )
        return bytes(blob)

    def _decompress_chunk(self, body_id: int, chunk_idx: int) -> bytes:
        """Decompress a single chunk (v2)."""
        key = (body_id, chunk_idx)
        # Lock-free fast path: chunk hits dominate the eval hot path
        # (decade-scale chunks vs day-scale queries), and dict.get is
        # atomic under the GIL — same argument as _eval_cache.
        cached = self._chunk_cache.get(key)
        if cached is not None:
            return cached
        # Double-checked lock: decompression is ~0.05 ms per chunk, so
        # holding the lock through the work is adequate and concurrent
        # requests for the same chunk never redo the zstd work.
        with self._decomp_lock:
            cached = self._chunk_cache.get(key)
            if cached is not None:
                return cached
            return self._decompress_chunk_uncached(body_id, chunk_idx)

    def _decompress_chunk_uncached(self, body_id: int, chunk_idx: int) -> bytes:
        """Decompress a single chunk (v2) and store it in the chunk cache."""
        key = (body_id, chunk_idx)
        body = self._bodies[body_id]
        chunk = self._chunk_index[body_id][chunk_idx]

        # The stored uncompressed size is fully derivable; validating it
        # caps the decompression allocation for corrupted/crafted files.
        expected = chunk.segment_count * (body.degree + 1) * body.components * 8
        if chunk.uncompressed_size != expected:
            raise LEBCorruptionError(
                f"Corrupted LEB2 chunk (body {body_id}, chunk {chunk_idx}): "
                f"uncompressed_size {chunk.uncompressed_size} != expected {expected}"
            )

        compressed = self._read_blob(
            chunk.blob_offset,
            chunk.compressed_size,
            f"chunk (body {body_id}, chunk {chunk_idx})",
        )
        decompressed = decompress_body(
            compressed,
            chunk.uncompressed_size,
            chunk.segment_count,
            body.degree,
            body.components,
        )

        # Bounded cache: keep max 64 chunks in memory. Runs under
        # _decomp_lock (held by the caller).
        if len(self._chunk_cache) > 64:
            self._chunk_cache.clear()
        self._chunk_cache[key] = decompressed
        return decompressed

    def _decompress_body(self, body_id: int) -> None:
        """Decompress a body's full coefficients (v1 legacy)."""
        # Lock-free fast path (see _decompress_chunk)
        if body_id in self._cache:
            return
        with self._decomp_lock:
            if body_id not in self._cache:
                self._decompress_body_uncached(body_id)

    def _decompress_body_uncached(self, body_id: int) -> bytes:
        """Decompress a v1 body's coefficients and store them in the cache."""
        entry = self._bodies[body_id]
        expected = entry.segment_count * (entry.degree + 1) * entry.components * 8
        if entry.uncompressed_size != expected:
            raise LEBCorruptionError(
                f"Corrupted LEB2 body entry {body_id}: uncompressed_size "
                f"{entry.uncompressed_size} != expected {expected}"
            )
        compressed = self._read_blob(
            entry.data_offset, entry.compressed_size, f"body {body_id}"
        )
        data = decompress_body(
            compressed,
            entry.uncompressed_size,
            entry.segment_count,
            entry.degree,
            entry.components,
        )
        # Runs under _decomp_lock (held by the caller)
        self._cache[body_id] = data
        return data

    def __enter__(self) -> "LEB2Reader":
        return self

    def __exit__(self, *args) -> None:
        self.close()

    @property
    def path(self) -> str:
        return self._path

    @property
    def jd_range(self) -> Tuple[float, float]:
        return (self._header.jd_start, self._header.jd_end)

    def warm(self, jd_start: float, jd_end: float) -> None:
        """Pre-fault mmap pages for data covering [jd_start, jd_end].

        Calls ``madvise(MADV_WILLNEED)`` on the byte ranges of chunks (v2)
        or body blobs (v1) that overlap the requested JD range.  Ranges are
        page-aligned and merged to minimise syscalls.

        Args:
            jd_start: Start of the Julian Day range to pre-fault.
            jd_end: End of the Julian Day range to pre-fault.
        """
        ranges: list[tuple[int, int]] = []
        if self._chunked:
            for chunks in self._chunk_index.values():
                for chunk in chunks:
                    if chunk.jd_end < jd_start or chunk.jd_start > jd_end:
                        continue
                    ranges.append((chunk.blob_offset, chunk.compressed_size))
        else:
            # v1 non-chunked: the body data is a single compressed blob,
            # so we cannot target individual segments within it.  Warm the
            # entire blob if the body's JD range overlaps.
            for entry in self._bodies.values():
                if entry.jd_end < jd_start or entry.jd_start > jd_end:
                    continue
                ranges.append((entry.data_offset, entry.compressed_size))

        _madvise_ranges(self._mm, ranges)

    def cool(self) -> None:
        """Advise the kernel that mmap pages can be reclaimed.

        Idempotent and safe to call on closed readers.
        """
        if self._mm is not None:
            _madvise_dontneed(self._mm)

    def has_body(self, body_id: int) -> bool:
        return body_id in self._bodies

    def eval_body(
        self, body_id: int, jd: float
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """Evaluate a body's position and velocity at a given Julian Day."""
        # Check eval cache
        _cache_key = (body_id, jd)
        cached = self._eval_cache.get(_cache_key)
        if cached is not None:
            return cached

        if body_id not in self._bodies:
            raise KeyError(f"Body {body_id} not in LEB file")

        body = self._bodies[body_id]

        if jd < body.jd_start or jd > body.jd_end:
            raise ValueError(
                f"JD {jd} outside range [{body.jd_start}, {body.jd_end}] "
                f"for body {body_id}"
            )

        # Corrupted-header guard (zero interval / empty segment table)
        if body.interval_days <= 0.0 or body.segment_count <= 0:
            raise ValueError(
                f"Corrupted LEB2 body entry {body_id}: interval_days="
                f"{body.interval_days}, segment_count={body.segment_count}"
            )

        # Global segment index
        seg_idx = int((jd - body.jd_start) / body.interval_days)
        seg_idx = max(0, min(seg_idx, body.segment_count - 1))

        # Get coefficients — chunked (v2) or monolithic (v1)
        if self._chunked:
            chunks = self._chunk_index[body_id]
            if not chunks:
                raise ValueError(
                    f"Corrupted LEB2 file: empty chunk index for body {body_id}"
                )
            chunk_idx = self._find_chunk(body_id, jd)
            chunk = chunks[chunk_idx]
            # The jd binary search compares against stored (write-time
            # rounded) chunk boundaries; in the half-ulp window where they
            # disagree with the segment arithmetic, trust the segment
            # index — clamping would silently evaluate the neighbouring
            # segment's polynomial.
            if not (
                chunk.segment_start
                <= seg_idx
                < chunk.segment_start + chunk.segment_count
            ):
                for step in (-1, 1):
                    alt = chunk_idx + step
                    if 0 <= alt < len(chunks):
                        cand = chunks[alt]
                        if (
                            cand.segment_start
                            <= seg_idx
                            < cand.segment_start + cand.segment_count
                        ):
                            chunk_idx, chunk = alt, cand
                            break
                else:
                    raise ValueError(
                        f"Corrupted LEB2 chunk index for body {body_id}: "
                        f"segment {seg_idx} not covered by any chunk near "
                        f"jd {jd}"
                    )
            chunk_data = self._decompress_chunk(body_id, chunk_idx)
            local_seg = seg_idx - chunk.segment_start
            deg1 = body.degree + 1
            n_coeffs = body.components * deg1
            byte_offset = local_seg * n_coeffs * 8
            coeffs = struct.unpack_from(f"<{n_coeffs}d", chunk_data, byte_offset)
        else:
            # v1: decompress full body on first access
            if body_id not in self._cache:
                self._decompress_body(body_id)
            deg1 = body.degree + 1
            n_coeffs = body.components * deg1
            byte_offset = seg_idx * n_coeffs * 8
            coeffs = struct.unpack_from(f"<{n_coeffs}d", self._cache[body_id], byte_offset)

        # Compute tau
        seg_start = body.jd_start + seg_idx * body.interval_days
        seg_mid = seg_start + 0.5 * body.interval_days
        tau = 2.0 * (jd - seg_mid) / body.interval_days
        if tau > 1.0:
            tau = 1.0
        elif tau < -1.0:
            tau = -1.0

        # Evaluate each component via Clenshaw
        pos = []
        vel = []
        scale = 2.0 / body.interval_days

        for c in range(body.components):
            comp_coeffs = coeffs[c * deg1: (c + 1) * deg1]
            val, deriv = _clenshaw_with_derivative(comp_coeffs, tau)
            pos.append(val)
            vel.append(deriv * scale)

        # Wrap longitude for ecliptic-frame bodies
        if body.coord_type in (COORD_ECLIPTIC, COORD_HELIO_ECL, COORD_GEO_ECLIPTIC):
            pos[0] = pos[0] % 360.0

        result = tuple(pos), tuple(vel)  # type: ignore[return-value]

        # Store in eval cache (bounded)
        if len(self._eval_cache) > 256:
            self._eval_cache.clear()
        self._eval_cache[_cache_key] = result  # type: ignore[assignment]

        return result  # type: ignore[return-value]

    def has_nutation(self) -> bool:
        """Return True if this LEB file contains usable nutation data.

        An empty (--skip-aux) nutation section has segment_count == 0;
        treating it as present would read the adjacent section's bytes as
        coefficients.
        """
        return self._nutation is not None and self._nutation.segment_count > 0

    def eval_nutation(self, jd_tt: float) -> Tuple[float, float]:
        """Evaluate nutation angles. Same as LEBReader.eval_nutation()."""
        if self._nutation is None or self._nutation.segment_count <= 0:
            raise ValueError("No nutation data in this LEB file")

        nut = self._nutation
        if nut.interval_days <= 0.0:
            raise ValueError("Corrupted LEB2 nutation header: zero interval")
        if jd_tt < nut.jd_start or jd_tt > nut.jd_end:
            raise ValueError(
                f"JD {jd_tt} outside nutation range [{nut.jd_start}, {nut.jd_end}]"
            )

        seg_idx = int((jd_tt - nut.jd_start) / nut.interval_days)
        seg_idx = max(0, min(seg_idx, nut.segment_count - 1))

        seg_start = nut.jd_start + seg_idx * nut.interval_days
        seg_mid = seg_start + 0.5 * nut.interval_days
        tau = 2.0 * (jd_tt - seg_mid) / nut.interval_days
        if tau > 1.0:
            tau = 1.0
        elif tau < -1.0:
            tau = -1.0

        deg1 = nut.degree + 1
        n_coeffs = nut.components * deg1
        seg_size = n_coeffs * 8
        byte_offset = self._nutation_data_offset + seg_idx * seg_size
        coeffs = struct.unpack_from(f"<{n_coeffs}d", self._mm, byte_offset)

        dpsi = _clenshaw(coeffs[0:deg1], tau)
        deps = _clenshaw(coeffs[deg1: 2 * deg1], tau)
        return dpsi, deps

    def delta_t(self, jd: float) -> float:
        """Get Delta-T. Same as LEBReader.delta_t()."""
        if not self._delta_t_jds:
            raise ValueError("No Delta-T data in this LEB file")

        jds = self._delta_t_jds
        vals = self._delta_t_vals
        n = len(jds)

        if jd <= jds[0]:
            return vals[0]
        if jd >= jds[-1]:
            return vals[-1]

        idx = bisect_right(jds, jd) - 1
        idx = max(0, min(idx, n - 2))

        span = jds[idx + 1] - jds[idx]
        if span == 0.0:  # pragma: no cover - sorted Delta-T table has no duplicate adjacent JDs
            return vals[idx]
        t = (jd - jds[idx]) / span
        return vals[idx] + t * (vals[idx + 1] - vals[idx])

    def get_star(self, star_id: int) -> StarEntry:
        """Look up a fixed star. Same as LEBReader.get_star()."""
        if star_id not in self._stars:
            raise KeyError(f"Star {star_id} not in LEB catalog")
        return self._stars[star_id]

    def close(self) -> None:
        """Close the memory-mapped file and release resources."""
        self._eval_cache.clear()
        self._cache.clear()
        self._chunk_cache.clear()
        self._chunk_index.clear()
        if self._mm is not None:
            try:
                self._mm.close()
            except (OSError, ValueError, KeyError):
                pass
            self._mm = None  # type: ignore[assignment]
        if self._file is not None:
            try:
                self._file.close()
            except (OSError, ValueError, KeyError):
                pass
            self._file = None  # type: ignore[assignment]
