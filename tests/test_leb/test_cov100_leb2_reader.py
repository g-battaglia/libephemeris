# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for ``libephemeris.leb2_reader``.

These tests drive the remaining uncovered lines and branch arcs in
:class:`~libephemeris.leb2_reader.LEB2Reader`: the missing-file guard,
parse-time error normalization (magic / version / compression scheme),
the optional-section presence branches, both the v2 chunked and v1
monolithic coefficient paths, the chunk neighbour search, the eval and
chunk cache eviction branches, the corrupted-size guards, the nutation
and Delta-T error guards plus their boundary returns, ``get_star``, and
the exception-swallowing arms of ``close``.

Every fixture is a small synthetic ``.leb2`` buffer built in ``tmp_path``
(no network, no external data files).  Genuine zstd-compressed blobs are
produced via :func:`libephemeris.leb_compression.compress_body` so the
decompression pipeline runs end-to-end.  Edge branches that cannot be
reached with a well-formed file are exercised by lightly mutating the
per-reader instance state only; the live module is never reloaded.
"""

from __future__ import annotations

import struct
from typing import Any

import numpy as np
import pytest

from libephemeris.leb_compression import compress_body
from libephemeris.leb_format import (
    CHUNK_ENTRY_SIZE,
    CHUNK_INDEX_HEADER_SIZE,
    COMPRESSED_BODY_ENTRY_SIZE,
    COMPRESSION_NONE,
    COMPRESSION_ZSTD_TRUNC_SHUFFLE,
    COORD_ECLIPTIC,
    COORD_ICRS_BARY,
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
    FileHeader,
    NutationHeader,
    SectionEntry,
    StarEntry,
    write_chunk_entry,
    write_chunk_index_header,
    write_compressed_body_entry,
    write_header,
    write_nutation_header,
    write_section_dir,
    write_star_entry,
)
from libephemeris.leb2_reader import LEB2Reader

JD0 = 2451545.0
SUN = 0
DEGREE = 2
COMPONENTS = 3
INTERVAL = 10.0


# =============================================================================
# RAW COEFFICIENT + COMPRESSED-BLOB HELPERS
# =============================================================================


def _raw_segments(n_seg: int, degree: int = DEGREE, components: int = COMPONENTS):
    """Build raw segment-major coefficient bytes for ``n_seg`` segments.

    Args:
        n_seg: Number of Chebyshev segments.
        degree: Polynomial degree.
        components: Components per segment.

    Returns:
        Bytes in segment-major float64 layout.
    """
    deg1 = degree + 1
    arr = np.zeros((n_seg, components, deg1), dtype=np.float64)
    # Non-trivial constant + linear terms so Clenshaw returns finite values.
    for s in range(n_seg):
        for c in range(components):
            arr[s, c, 0] = float(s + c + 1)
            if deg1 > 1:
                arr[s, c, 1] = 0.5
    return arr.tobytes()


def _compress(raw: bytes, n_seg: int, degree: int = DEGREE, components: int = COMPONENTS):
    """Compress raw coefficients with the zstd-truncate-shuffle pipeline."""
    deg1 = degree + 1
    bits = [52] * deg1  # keep all mantissa bits — lossless round trip
    return compress_body(raw, n_seg, degree, components, bits)


# =============================================================================
# SYNTHETIC FILE BUILDERS
# =============================================================================


def _build_v2_file(
    path,
    *,
    n_seg: int = 30,
    chunk_segments: int = 10,
    interval: float = INTERVAL,
    coord_type: int = COORD_ICRS_BARY,
    with_nutation: bool = False,
    with_delta_t: bool = False,
    with_stars: bool = False,
    flags: int = COMPRESSION_ZSTD_TRUNC_SHUFFLE,
    version: int = LEB2_VERSION,
    magic: bytes = LEB2_MAGIC,
) -> str:
    """Build a valid v2 chunked LEB2 file with one body (SUN).

    The body is split into chunks of ``chunk_segments`` segments each.
    Optional nutation / Delta-T / star sections may be appended so the
    reader exposes the corresponding mutable state.

    Args:
        path: Destination path object.
        n_seg: Total segment count for the body.
        chunk_segments: Segments per chunk.
        interval: Segment interval in days.
        coord_type: Body coordinate type (selects longitude-wrap branch).
        with_nutation: Append a nutation section.
        with_delta_t: Append a Delta-T section.
        with_stars: Append a star catalog section.
        flags: Header compression-scheme flag.
        version: Header version field.
        magic: Header magic bytes.

    Returns:
        The string path to the written file.
    """
    deg1 = DEGREE + 1
    n_coeffs = COMPONENTS * deg1
    seg_bytes = n_coeffs * 8

    # Build the per-chunk compressed blobs.
    raw = _raw_segments(n_seg)
    chunk_specs = []  # (segment_start, segment_count, compressed_blob, uncompressed)
    seg_start = 0
    while seg_start < n_seg:
        cnt = min(chunk_segments, n_seg - seg_start)
        chunk_raw = raw[seg_start * seg_bytes : (seg_start + cnt) * seg_bytes]
        blob = _compress(chunk_raw, cnt)
        chunk_specs.append((seg_start, cnt, blob, cnt * seg_bytes))
        seg_start += cnt

    n_sections = 1 + int(with_nutation) + int(with_delta_t) + int(with_stars)
    sec_dir_off = HEADER_SIZE
    body_idx_off = sec_dir_off + n_sections * SECTION_DIR_SIZE
    body_data_off = body_idx_off + COMPRESSED_BODY_ENTRY_SIZE

    # Body data = chunk index header + chunk entries + concatenated blobs.
    chunk_index_off = body_data_off
    chunk_entries_off = chunk_index_off + CHUNK_INDEX_HEADER_SIZE
    blobs_off = chunk_entries_off + len(chunk_specs) * CHUNK_ENTRY_SIZE

    # Lay out blob offsets.
    blob_layout = []  # (blob_offset, blob)
    cursor = blobs_off
    for (_ss, _cnt, blob, _u) in chunk_specs:
        blob_layout.append((cursor, blob))
        cursor += len(blob)

    body_data_size = cursor - body_data_off

    sections: list[tuple[int, int, int]] = []
    sections.append((SECTION_BODY_INDEX, body_idx_off, COMPRESSED_BODY_ENTRY_SIZE))

    nut_off = nut_data_off = 0
    nut_nseg = nut_deg = nut_comps = 0
    nut_interval = 20.0
    if with_nutation:
        nut_deg, nut_comps, nut_nseg = 2, 2, 2
        nut_data_size = nut_nseg * nut_comps * (nut_deg + 1) * 8
        nut_off = cursor
        nut_data_off = nut_off + NUTATION_HEADER_SIZE
        nut_size = NUTATION_HEADER_SIZE + nut_data_size
        sections.append((SECTION_NUTATION, nut_off, nut_size))
        cursor += nut_size

    dt_off = 0
    dt_entries: list[tuple[float, float]] = []
    if with_delta_t:
        dt_entries = [
            (JD0, 0.0008),
            (JD0 + 30.0, 0.0008),  # equal value -> span!=0 still, see span==0 test
            (JD0 + 60.0, 0.0010),
        ]
        dt_off = cursor
        dt_size = DELTA_T_HEADER_SIZE + len(dt_entries) * DELTA_T_ENTRY_SIZE
        sections.append((SECTION_DELTA_T, dt_off, dt_size))
        cursor += dt_size

    star_off = 0
    star_ids: list[int] = []
    if with_stars:
        star_ids = [1, 2]
        star_off = cursor
        star_size = len(star_ids) * STAR_ENTRY_SIZE
        sections.append((SECTION_STARS, star_off, star_size))
        cursor += star_size

    buf = bytearray(cursor)
    write_header(
        buf,
        FileHeader(
            magic=magic,
            version=version,
            section_count=n_sections,
            body_count=1,
            jd_start=JD0,
            jd_end=JD0 + n_seg * interval,
            generation_epoch=0.0,
            flags=flags,
        ),
    )
    for i, (sid, off, size) in enumerate(sections):
        write_section_dir(
            buf, sec_dir_off + i * SECTION_DIR_SIZE, SectionEntry(sid, off, size)
        )

    write_compressed_body_entry(
        buf,
        body_idx_off,
        CompressedBodyEntry(
            body_id=SUN,
            coord_type=coord_type,
            segment_count=n_seg,
            jd_start=JD0,
            jd_end=JD0 + n_seg * interval,
            interval_days=interval,
            degree=DEGREE,
            components=COMPONENTS,
            data_offset=body_data_off,
            compressed_size=body_data_size,
            uncompressed_size=n_seg * seg_bytes,
        ),
    )

    # Chunk index header + entries.
    write_chunk_index_header(buf, chunk_index_off, len(chunk_specs), interval)
    for i, (ss, cnt, blob, unc) in enumerate(chunk_specs):
        boff, _b = blob_layout[i]
        write_chunk_entry(
            buf,
            chunk_entries_off + i * CHUNK_ENTRY_SIZE,
            ChunkEntry(
                jd_start=JD0 + ss * interval,
                jd_end=JD0 + (ss + cnt) * interval,
                blob_offset=boff,
                compressed_size=len(blob),
                uncompressed_size=unc,
                segment_start=ss,
                segment_count=cnt,
            ),
        )
    for (boff, blob) in blob_layout:
        buf[boff : boff + len(blob)] = blob

    if with_nutation:
        write_nutation_header(
            buf,
            nut_off,
            NutationHeader(
                jd_start=JD0,
                jd_end=JD0 + nut_nseg * nut_interval,
                interval_days=nut_interval,
                degree=nut_deg,
                components=nut_comps,
                segment_count=nut_nseg,
                reserved=0,
            ),
        )
        stride = nut_comps * (nut_deg + 1) * 8
        for s in range(nut_nseg):
            for ci in range(nut_comps * (nut_deg + 1)):
                struct.pack_into(
                    "<d", buf, nut_data_off + s * stride + ci * 8, 1e-5 * (ci + 1)
                )

    if with_delta_t:
        struct.pack_into(DELTA_T_HEADER_FMT, buf, dt_off, len(dt_entries), 0)
        data_off = dt_off + DELTA_T_HEADER_SIZE
        for i, (jd, dt) in enumerate(dt_entries):
            struct.pack_into(
                DELTA_T_ENTRY_FMT, buf, data_off + i * DELTA_T_ENTRY_SIZE, jd, dt
            )

    if with_stars:
        for i, sid in enumerate(star_ids):
            write_star_entry(
                buf,
                star_off + i * STAR_ENTRY_SIZE,
                StarEntry(
                    star_id=sid,
                    ra_j2000=10.0 * sid,
                    dec_j2000=5.0 * sid,
                    pm_ra=0.0,
                    pm_dec=0.0,
                    parallax=0.0,
                    rv=0.0,
                    magnitude=1.0,
                ),
            )

    path.write_bytes(bytes(buf))
    return str(path)


def _build_v1_file(
    path,
    *,
    n_seg: int = 5,
    interval: float = INTERVAL,
    coord_type: int = COORD_ICRS_BARY,
) -> str:
    """Build a valid v1 monolithic LEB2 file with one body (SUN).

    Args:
        path: Destination path object.
        n_seg: Segment count for the body.
        interval: Segment interval in days.
        coord_type: Body coordinate type.

    Returns:
        The string path to the written file.
    """
    deg1 = DEGREE + 1
    n_coeffs = COMPONENTS * deg1
    seg_bytes = n_coeffs * 8

    raw = _raw_segments(n_seg)
    blob = _compress(raw, n_seg)

    n_sections = 1
    sec_dir_off = HEADER_SIZE
    body_idx_off = sec_dir_off + n_sections * SECTION_DIR_SIZE
    body_data_off = body_idx_off + COMPRESSED_BODY_ENTRY_SIZE
    cursor = body_data_off + len(blob)

    buf = bytearray(cursor)
    write_header(
        buf,
        FileHeader(
            magic=LEB2_MAGIC,
            version=LEB2_VERSION_V1,
            section_count=n_sections,
            body_count=1,
            jd_start=JD0,
            jd_end=JD0 + n_seg * interval,
            generation_epoch=0.0,
            flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
        ),
    )
    write_section_dir(
        buf,
        sec_dir_off,
        SectionEntry(
            SECTION_BODY_INDEX, body_idx_off, COMPRESSED_BODY_ENTRY_SIZE
        ),
    )
    write_compressed_body_entry(
        buf,
        body_idx_off,
        CompressedBodyEntry(
            body_id=SUN,
            coord_type=coord_type,
            segment_count=n_seg,
            jd_start=JD0,
            jd_end=JD0 + n_seg * interval,
            interval_days=interval,
            degree=DEGREE,
            components=COMPONENTS,
            data_offset=body_data_off,
            compressed_size=len(blob),
            uncompressed_size=n_seg * seg_bytes,
        ),
    )
    buf[body_data_off : body_data_off + len(blob)] = blob
    path.write_bytes(bytes(buf))
    return str(path)


class _JDStartTrap:
    """Proxy over a body/nutation entry that flips ``jd_start`` per access.

    The first ``jd_start`` read (the range check) returns a *low* value so the
    in-range guard passes; subsequent reads (the segment-origin / tau math)
    return a *high* value so ``seg_mid`` sits above ``jd`` and the computed
    tau falls below -1.0, forcing the ``tau < -1.0`` clamp branch.  Every
    other attribute delegates to the wrapped entry unchanged.
    """

    def __init__(self, base: Any, delta: float = 5000.0) -> None:
        object.__setattr__(self, "_base", base)
        object.__setattr__(self, "_delta", delta)
        object.__setattr__(self, "_n", 0)

    @property
    def jd_start(self) -> float:  # type: ignore[override]
        base = object.__getattribute__(self, "_base")
        delta = object.__getattribute__(self, "_delta")
        n = object.__getattribute__(self, "_n") + 1
        object.__setattr__(self, "_n", n)
        return base.jd_start - delta if n == 1 else base.jd_start + delta

    def __getattr__(self, name: str) -> Any:
        return getattr(object.__getattribute__(self, "_base"), name)


# =============================================================================
# __init__ guards + parse error normalization
# =============================================================================


class TestInitAndParseErrors:
    """Cover the missing-file guard and parse-time error normalization."""

    @pytest.mark.unit
    def test_missing_file_raises(self, tmp_path):
        """A non-existent path raises FileNotFoundError (line 78)."""
        with pytest.raises(FileNotFoundError, match="LEB file not found"):
            LEB2Reader(str(tmp_path / "does_not_exist.leb2"))

    @pytest.mark.unit
    def test_truncated_file_normalized_to_value_error(self, tmp_path):
        """A header claiming sections but lacking their bytes -> ValueError.

        Drives the ``except (struct.error, IndexError, KeyError)`` arm that
        normalizes truncated reads into a single ValueError (lines 94-98).
        """
        buf = bytearray(HEADER_SIZE)
        write_header(
            buf,
            FileHeader(
                magic=LEB2_MAGIC,
                version=LEB2_VERSION,
                section_count=5,  # claims 5 sections; bytes absent -> struct.error
                body_count=1,
                jd_start=JD0,
                jd_end=JD0 + 100.0,
                generation_epoch=0.0,
                flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
            ),
        )
        path = tmp_path / "trunc.leb2"
        path.write_bytes(bytes(buf))
        with pytest.raises(ValueError, match="Corrupted LEB2 file"):
            LEB2Reader(str(path))

    @pytest.mark.unit
    def test_value_error_reraised_and_closed(self, tmp_path):
        """An invalid magic -> ValueError re-raised by the outer arm.

        Covers the ``except (OSError, ValueError)`` re-raise path
        (lines 99-101) plus the invalid-magic guard (line 106).
        """
        buf = bytearray(HEADER_SIZE)
        write_header(
            buf,
            FileHeader(
                magic=b"XXXX",
                version=LEB2_VERSION,
                section_count=0,
                body_count=0,
                jd_start=JD0,
                jd_end=JD0 + 100.0,
                generation_epoch=0.0,
                flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
            ),
        )
        path = tmp_path / "badmagic.leb2"
        path.write_bytes(bytes(buf))
        with pytest.raises(ValueError, match="Invalid LEB2 magic"):
            LEB2Reader(str(path))

    @pytest.mark.unit
    def test_unsupported_version_raises(self, tmp_path):
        """A magic-valid file with an unknown version -> ValueError (line 110)."""
        buf = bytearray(HEADER_SIZE)
        write_header(
            buf,
            FileHeader(
                magic=LEB2_MAGIC,
                version=99,
                section_count=0,
                body_count=0,
                jd_start=JD0,
                jd_end=JD0 + 100.0,
                generation_epoch=0.0,
                flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
            ),
        )
        path = tmp_path / "badver.leb2"
        path.write_bytes(bytes(buf))
        with pytest.raises(ValueError, match="Unsupported LEB2 version"):
            LEB2Reader(str(path))

    @pytest.mark.unit
    def test_unsupported_compression_scheme_raises(self, tmp_path):
        """An unknown compression flag -> ValueError (line 122)."""
        buf = bytearray(HEADER_SIZE)
        write_header(
            buf,
            FileHeader(
                magic=LEB2_MAGIC,
                version=LEB2_VERSION,
                section_count=0,
                body_count=0,
                jd_start=JD0,
                jd_end=JD0 + 100.0,
                generation_epoch=0.0,
                flags=COMPRESSION_NONE,
            ),
        )
        path = tmp_path / "badcomp.leb2"
        path.write_bytes(bytes(buf))
        with pytest.raises(ValueError, match="Unsupported LEB2 compression scheme"):
            LEB2Reader(str(path))


# =============================================================================
# Optional-section presence branches in _parse
# =============================================================================


class TestParseOptionalSections:
    """Cover the section presence / absence branch arcs in ``_parse``."""

    @pytest.mark.unit
    def test_all_optional_sections_present(self, tmp_path):
        """A file with body+nutation+delta-t+stars takes every *true* arc.

        Covers 135->143 (body index present), the chunked-index loop,
        150->156 (nutation present), 158->169 (delta-t present),
        170->exit (stars present).
        """
        path = _build_v2_file(
            tmp_path / "full.leb2",
            with_nutation=True,
            with_delta_t=True,
            with_stars=True,
        )
        reader = LEB2Reader(path)
        try:
            assert reader.has_body(SUN)
            assert reader._nutation is not None
            assert reader._delta_t_jds  # populated
            assert reader._stars  # populated
        finally:
            reader.close()

    @pytest.mark.unit
    def test_no_optional_sections(self, tmp_path):
        """A v1 body-only file takes every *false* optional-section arc.

        Covers 143->148 (not chunked -> skip chunk index), the absent
        nutation / delta-t / stars false arcs, and 170->exit fall-through.
        """
        path = _build_v1_file(tmp_path / "bodyonly.leb2")
        reader = LEB2Reader(path)
        try:
            assert reader.has_body(SUN)
            assert reader._nutation is None
            assert reader._delta_t_jds == []
            assert reader._stars == {}
        finally:
            reader.close()

    @pytest.mark.unit
    def test_no_body_index_section(self, tmp_path):
        """A v2 file with no body-index section skips body parsing (135->143).

        The section directory is empty, so ``SECTION_BODY_INDEX not in
        self._sections`` jumps straight to the chunked guard (line 143)
        with an empty ``_bodies`` dict.
        """
        buf = bytearray(HEADER_SIZE)
        write_header(
            buf,
            FileHeader(
                magic=LEB2_MAGIC,
                version=LEB2_VERSION,  # chunked
                section_count=0,
                body_count=0,
                jd_start=JD0,
                jd_end=JD0 + 100.0,
                generation_epoch=0.0,
                flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
            ),
        )
        path = tmp_path / "nobody.leb2"
        path.write_bytes(bytes(buf))
        reader = LEB2Reader(str(path))
        try:
            assert reader._bodies == {}
            assert reader._chunk_index == {}
            assert reader._chunked is True
        finally:
            reader.close()


# =============================================================================
# eval_body — v2 chunked path, cache, guards, neighbour search
# =============================================================================


class TestEvalBodyChunked:
    """Cover the v2 chunked eval path and its guards/caches."""

    @pytest.mark.unit
    def test_eval_and_cache_hit(self, tmp_path):
        """A repeated eval hits the eval cache (line 321) and chunk cache (208)."""
        path = _build_v2_file(tmp_path / "ev.leb2")
        reader = LEB2Reader(path)
        try:
            jd = JD0 + 5.0
            r1 = reader.eval_body(SUN, jd)
            r2 = reader.eval_body(SUN, jd)  # eval cache hit -> line 321
            assert r1 == r2
            # Force chunk cache hit (line 208): a different jd in the same
            # chunk re-reads the cached decompressed chunk.
            r3 = reader.eval_body(SUN, jd + 1.0)
            assert len(r3[0]) == 3
        finally:
            reader.close()

    @pytest.mark.unit
    def test_body_not_found(self, tmp_path):
        """An unknown body raises KeyError (line 324)."""
        path = _build_v2_file(tmp_path / "nf.leb2")
        reader = LEB2Reader(path)
        try:
            with pytest.raises(KeyError, match="not in LEB file"):
                reader.eval_body(999, JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_jd_out_of_range(self, tmp_path):
        """A JD outside the body range raises ValueError (line 329)."""
        path = _build_v2_file(tmp_path / "oor.leb2")
        reader = LEB2Reader(path)
        try:
            with pytest.raises(ValueError, match="outside range"):
                reader.eval_body(SUN, JD0 - 100.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_corrupted_interval_guard(self, tmp_path):
        """Zero interval -> corrupted-header ValueError (line 336)."""
        path = _build_v2_file(tmp_path / "corr.leb2")
        reader = LEB2Reader(path)
        try:
            reader._bodies[SUN].interval_days = 0.0
            with pytest.raises(ValueError, match="Corrupted LEB2 body entry"):
                reader.eval_body(SUN, JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_empty_chunk_index_guard(self, tmp_path):
        """An empty chunk index -> ValueError (line 349)."""
        path = _build_v2_file(tmp_path / "emptychunks.leb2")
        reader = LEB2Reader(path)
        try:
            reader._chunk_index[SUN] = []
            with pytest.raises(ValueError, match="empty chunk index"):
                reader.eval_body(SUN, JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_neighbour_chunk_search_finds_segment(self, tmp_path):
        """A mismatched chunk boundary triggers the neighbour search (364-374).

        The binary search lands on a chunk whose segment range does not
        contain ``seg_idx``; the +/-1 neighbour scan recovers the correct
        chunk instead of raising.
        """
        path = _build_v2_file(
            tmp_path / "neigh.leb2", n_seg=30, chunk_segments=10
        )
        reader = LEB2Reader(path)
        try:
            chunks = reader._chunk_index[SUN]
            assert len(chunks) >= 2
            # Corrupt the first chunk's jd_end so the binary search for a JD
            # in chunk 1 instead settles on chunk 0; the segment index still
            # belongs to chunk 1, so the +1 neighbour scan recovers it.
            target_jd = JD0 + 15.0 * INTERVAL / INTERVAL  # within chunk 1
            target_jd = JD0 + 12.0 * INTERVAL  # seg 12 -> chunk index 1
            # Make chunk 0 claim it extends past the target so _find_chunk
            # returns 0, but its segment range (0..9) excludes seg 12.
            chunks[0].jd_end = chunks[1].jd_end + 1.0
            pos, vel = reader.eval_body(SUN, target_jd)
            assert len(pos) == 3 and len(vel) == 3
        finally:
            reader.close()

    @pytest.mark.unit
    def test_neighbour_chunk_search_not_found(self, tmp_path):
        """No covering neighbour chunk -> corrupted-index ValueError (line 376)."""
        path = _build_v2_file(
            tmp_path / "notfound.leb2", n_seg=30, chunk_segments=10
        )
        reader = LEB2Reader(path)
        try:
            chunks = reader._chunk_index[SUN]
            # Break every chunk's segment range so the seg_idx falls in no
            # chunk; the binary search lands somewhere and neither neighbour
            # covers it.
            for ch in chunks:
                ch.segment_start = 10_000
                ch.segment_count = 1
            with pytest.raises(ValueError, match="not covered by any chunk"):
                reader.eval_body(SUN, JD0 + 12.0 * INTERVAL)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_tau_clamped_high(self, tmp_path):
        """Extending jd_end past real coverage forces tau > 1.0 (line 401)."""
        path = _build_v2_file(tmp_path / "tauhi.leb2", n_seg=10, chunk_segments=10)
        reader = LEB2Reader(path)
        try:
            body = reader._bodies[SUN]
            # Declare the body range far past the last segment so a JD near
            # jd_end clamps onto the final segment with tau > 1.
            body.jd_end = (
                body.jd_start + body.segment_count * body.interval_days + 1000.0
            )
            pos, vel = reader.eval_body(SUN, body.jd_end - 1.0)
            assert len(pos) == 3 and len(vel) == 3
        finally:
            reader.close()

    @pytest.mark.unit
    def test_tau_clamped_low(self, tmp_path):
        """A jd_start trap forces tau < -1.0 in eval_body (line 403)."""
        path = _build_v2_file(tmp_path / "taulo.leb2", n_seg=10, chunk_segments=10)
        reader = LEB2Reader(path)
        try:
            real = reader._bodies[SUN]
            reader._bodies[SUN] = _JDStartTrap(real)
            # Range check sees a low jd_start (passes); the segment-origin
            # math sees a high jd_start (seg_mid above jd) -> tau < -1.0.
            pos, vel = reader.eval_body(SUN, real.jd_start)
            assert len(pos) == 3 and len(vel) == 3
        finally:
            reader._bodies[SUN] = real
            reader.close()

    @pytest.mark.unit
    def test_eval_cache_eviction(self, tmp_path):
        """More than 256 distinct evals clears the eval cache (line 424)."""
        path = _build_v2_file(tmp_path / "evict.leb2", n_seg=30, chunk_segments=10)
        reader = LEB2Reader(path)
        try:
            span = reader._bodies[SUN].jd_end - reader._bodies[SUN].jd_start
            for i in range(300):
                jd = JD0 + (i / 300.0) * span * 0.999
                reader.eval_body(SUN, jd)
            # After eviction the cache is far below the count of queries.
            assert len(reader._eval_cache) <= 257
        finally:
            reader.close()

    @pytest.mark.unit
    def test_longitude_wrap_for_ecliptic(self, tmp_path):
        """Ecliptic-frame bodies wrap longitude into [0, 360) (line 418)."""
        path = _build_v2_file(
            tmp_path / "ecl.leb2", coord_type=COORD_ECLIPTIC
        )
        reader = LEB2Reader(path)
        try:
            pos, _vel = reader.eval_body(SUN, JD0 + 5.0)
            assert 0.0 <= pos[0] < 360.0
        finally:
            reader.close()


# =============================================================================
# eval_body — v1 monolithic path
# =============================================================================


class TestEvalBodyMonolithic:
    """Cover the v1 monolithic decompression path."""

    @pytest.mark.unit
    def test_v1_eval_decompresses_full_body(self, tmp_path):
        """A v1 file decompresses the whole body on first access (389-394)."""
        path = _build_v1_file(tmp_path / "v1.leb2", n_seg=5)
        reader = LEB2Reader(path)
        try:
            assert SUN not in reader._cache
            pos, vel = reader.eval_body(SUN, JD0 + 5.0)
            assert SUN in reader._cache  # full-body cache populated
            assert len(pos) == 3 and len(vel) == 3
            # Second access reuses the cache.
            reader.eval_body(SUN, JD0 + 15.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_v1_corrupted_size_guard(self, tmp_path):
        """A v1 body with a wrong uncompressed_size -> ValueError (241-247)."""
        path = _build_v1_file(tmp_path / "v1corr.leb2", n_seg=5)
        reader = LEB2Reader(path)
        try:
            reader._bodies[SUN].uncompressed_size += 8  # mismatch expected
            with pytest.raises(ValueError, match="Corrupted LEB2 body entry"):
                reader.eval_body(SUN, JD0 + 5.0)
        finally:
            reader.close()


# =============================================================================
# _decompress_chunk guards
# =============================================================================


class TestDecompressChunk:
    """Cover the chunk decompression guards and cache eviction."""

    @pytest.mark.unit
    def test_chunk_corrupted_size_guard(self, tmp_path):
        """A chunk with a wrong uncompressed_size -> ValueError (line 217)."""
        path = _build_v2_file(tmp_path / "chunkcorr.leb2")
        reader = LEB2Reader(path)
        try:
            reader._chunk_index[SUN][0].uncompressed_size += 8
            with pytest.raises(ValueError, match="Corrupted LEB2 chunk"):
                reader._decompress_chunk(SUN, 0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_chunk_cache_eviction(self, tmp_path):
        """More than 64 cached chunks clears the chunk cache (line 235)."""
        path = _build_v2_file(tmp_path / "chunkevict.leb2")
        reader = LEB2Reader(path)
        try:
            # Seed >64 dummy chunk-cache entries, then trigger a real
            # decompress so the size-guard passes and the clear branch runs.
            for i in range(70):
                reader._chunk_cache[(SUN, 1000 + i)] = b""
            assert len(reader._chunk_cache) > 64
            reader._decompress_chunk(SUN, 0)
            # The clear ran; only the freshly decompressed chunk remains.
            assert len(reader._chunk_cache) == 1
        finally:
            reader.close()


# =============================================================================
# path property + warm v1 branch
# =============================================================================


class TestPathAndWarm:
    """Cover the path property and the v1 warm branch."""

    @pytest.mark.unit
    def test_path_property(self, tmp_path):
        """The path property returns the file path (line 267)."""
        path = _build_v2_file(tmp_path / "pp.leb2")
        reader = LEB2Reader(path)
        try:
            assert reader.path == path
        finally:
            reader.close()

    @pytest.mark.unit
    def test_warm_v1_branch(self, tmp_path):
        """warm() on a v1 (non-chunked) reader walks body blobs (295-298)."""
        path = _build_v1_file(tmp_path / "warmv1.leb2", n_seg=5)
        reader = LEB2Reader(path)
        try:
            body = reader._bodies[SUN]
            # Overlapping range -> append branch (296->298 false, 298 taken).
            reader.warm(body.jd_start, body.jd_end)
            # Non-overlapping range -> the continue arc (296 true).
            reader.warm(body.jd_end + 1000.0, body.jd_end + 2000.0)
        finally:
            reader.close()


# =============================================================================
# eval_nutation guards + tau clamps
# =============================================================================


class TestNutation:
    """Cover eval_nutation error guards and tau clamps."""

    @pytest.mark.unit
    def test_no_nutation_data(self, tmp_path):
        """A reader without nutation -> ValueError (line 441)."""
        path = _build_v1_file(tmp_path / "nonut.leb2")
        reader = LEB2Reader(path)
        try:
            with pytest.raises(ValueError, match="No nutation data"):
                reader.eval_nutation(JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_zero_interval_nutation(self, tmp_path):
        """Zero nutation interval -> corrupted-header ValueError (line 445)."""
        path = _build_v2_file(tmp_path / "nutzero.leb2", with_nutation=True)
        reader = LEB2Reader(path)
        try:
            reader._nutation.interval_days = 0.0
            with pytest.raises(ValueError, match="zero interval"):
                reader.eval_nutation(JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_jd_out_of_nutation_range(self, tmp_path):
        """A JD before the nutation start -> ValueError (line 447)."""
        path = _build_v2_file(tmp_path / "nutrange.leb2", with_nutation=True)
        reader = LEB2Reader(path)
        try:
            nut = reader._nutation
            with pytest.raises(ValueError, match="outside nutation range"):
                reader.eval_nutation(nut.jd_start - 10.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_nutation_tau_clamped_high(self, tmp_path):
        """Extending nut jd_end forces tau > 1.0 (line 458)."""
        path = _build_v2_file(tmp_path / "nuttauhi.leb2", with_nutation=True)
        reader = LEB2Reader(path)
        try:
            nut = reader._nutation
            nut.jd_end = nut.jd_start + nut.segment_count * nut.interval_days + 500.0
            dpsi, deps = reader.eval_nutation(nut.jd_end - 1.0)
            assert isinstance(dpsi, float) and isinstance(deps, float)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_nutation_tau_clamped_low(self, tmp_path):
        """A jd_start trap forces tau < -1.0 in eval_nutation (line 460)."""
        path = _build_v2_file(tmp_path / "nuttaulo.leb2", with_nutation=True)
        reader = LEB2Reader(path)
        try:
            real = reader._nutation
            reader._nutation = _JDStartTrap(real)
            dpsi, deps = reader.eval_nutation(real.jd_start)
            assert isinstance(dpsi, float) and isinstance(deps, float)
        finally:
            reader._nutation = real
            reader.close()


# =============================================================================
# delta_t guards + boundary returns
# =============================================================================


class TestDeltaT:
    """Cover the delta_t no-data guard, boundary returns, and span==0."""

    @pytest.mark.unit
    def test_no_delta_t_data(self, tmp_path):
        """A reader without a Delta-T table -> ValueError (line 475)."""
        path = _build_v1_file(tmp_path / "nodt.leb2")
        reader = LEB2Reader(path)
        try:
            with pytest.raises(ValueError, match="No Delta-T data"):
                reader.delta_t(JD0 + 5.0)
        finally:
            reader.close()

    @pytest.mark.unit
    def test_delta_t_below_first(self, tmp_path):
        """A JD at/below the first sample returns vals[0] (line 482)."""
        path = _build_v2_file(tmp_path / "dtlow.leb2", with_delta_t=True)
        reader = LEB2Reader(path)
        try:
            first_jd = reader._delta_t_jds[0]
            assert reader.delta_t(first_jd - 100.0) == reader._delta_t_vals[0]
        finally:
            reader.close()

    @pytest.mark.unit
    def test_delta_t_above_last(self, tmp_path):
        """A JD at/above the last sample returns vals[-1] (line 484)."""
        path = _build_v2_file(tmp_path / "dthigh.leb2", with_delta_t=True)
        reader = LEB2Reader(path)
        try:
            last_jd = reader._delta_t_jds[-1]
            assert reader.delta_t(last_jd + 100.0) == reader._delta_t_vals[-1]
        finally:
            reader.close()

    @pytest.mark.unit
    def test_delta_t_interpolates(self, tmp_path):
        """A mid-range JD interpolates between samples (line 493 fall-through)."""
        path = _build_v2_file(tmp_path / "dtinterp.leb2", with_delta_t=True)
        reader = LEB2Reader(path)
        try:
            jds = reader._delta_t_jds
            mid = (jds[-2] + jds[-1]) / 2.0
            val = reader.delta_t(mid)
            lo = min(reader._delta_t_vals[-2], reader._delta_t_vals[-1])
            hi = max(reader._delta_t_vals[-2], reader._delta_t_vals[-1])
            assert lo <= val <= hi
        finally:
            reader.close()


# =============================================================================
# get_star
# =============================================================================


class TestGetStar:
    """Cover get_star hit and miss branches."""

    @pytest.mark.unit
    def test_get_star_found(self, tmp_path):
        """A present star is returned (line 499)."""
        path = _build_v2_file(tmp_path / "starhit.leb2", with_stars=True)
        reader = LEB2Reader(path)
        try:
            star = reader.get_star(1)
            assert star.star_id == 1
        finally:
            reader.close()

    @pytest.mark.unit
    def test_get_star_missing(self, tmp_path):
        """A missing star raises KeyError (lines 497-498)."""
        path = _build_v2_file(tmp_path / "starmiss.leb2", with_stars=True)
        reader = LEB2Reader(path)
        try:
            with pytest.raises(KeyError, match="not in LEB catalog"):
                reader.get_star(99999)
        finally:
            reader.close()


# =============================================================================
# close() exception swallowing
# =============================================================================


class TestCloseSwallowsErrors:
    """Cover the close() exception-swallowing arms (507->513, 510-511, 516-517)."""

    @pytest.mark.unit
    def test_close_ignores_mm_and_file_errors(self, tmp_path):
        """close() swallows errors from both mmap.close and file.close."""
        path = _build_v2_file(tmp_path / "close.leb2")
        reader = LEB2Reader(path)

        class _BadMM:
            def close(self) -> None:
                raise OSError("mm close failed")

        class _BadFile:
            def close(self) -> None:
                raise ValueError("file close failed")

        # Release the real handles first to free OS resources, then swap in
        # raising stand-ins to exercise the except arms (510-511, 516-517).
        reader._mm.close()
        reader._file.close()
        reader._mm = _BadMM()  # type: ignore[assignment]
        reader._file = _BadFile()  # type: ignore[assignment]

        reader.close()  # must not raise
        assert reader._mm is None
        assert reader._file is None

    @pytest.mark.unit
    def test_close_with_none_handles(self, tmp_path):
        """close() with already-None handles takes the false guard arcs.

        Covers 507->513 (mm is None -> skip) and 513->exit (file is None).
        """
        path = _build_v2_file(tmp_path / "closenone.leb2")
        reader = LEB2Reader(path)
        reader.close()  # first close nulls both handles
        reader.close()  # second close: both None -> guard-false arcs
        assert reader._mm is None
        assert reader._file is None
