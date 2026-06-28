# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
LEB (LibEphemeris Binary) file format definitions.

This module defines the binary file format constants, dataclasses, and
serialization helpers used by both the generator (write) and reader (read).

The .leb format stores precomputed Chebyshev polynomial coefficients for
fast ephemeris evaluation, providing 10-20x speedup over the Skyfield pipeline.

File layout:
    0x0000       File Header (64 bytes)
    0x0040       Section Directory (N x 24 bytes)
    variable     Section 0: Body Index
    variable     Section 1: Chebyshev Coefficient Data
    variable     Section 2: Nutation Chebyshev Data
    variable     Section 3: Delta-T Sparse Table
    variable     Section 4: Star Catalog
    variable     Section 5: Orbital Elements (reserved)
"""

from __future__ import annotations

import os
import struct
from dataclasses import dataclass
from typing import Any

# =============================================================================
# FORMAT CONSTANTS
# =============================================================================

MAGIC = b"LEB1"
VERSION = 1

# Coordinate types
COORD_ICRS_BARY = 0  # ICRS barycentric planet center (x, y, z) in AU
COORD_ECLIPTIC = 1  # Ecliptic of date (lon, lat, dist) in deg/deg/AU
COORD_HELIO_ECL = 2  # Heliocentric ecliptic (lon, lat, dist) in deg/deg/AU
COORD_GEO_ECLIPTIC = 3  # Reserved — not used (retrograde cusps; see algorithms.md §16)
COORD_ICRS_BARY_SYSTEM = 4  # ICRS system barycenter (x, y, z) in AU — runtime COB

# Section IDs
SECTION_BODY_INDEX = 0
SECTION_CHEBYSHEV = 1
SECTION_NUTATION = 2
SECTION_DELTA_T = 3
SECTION_STARS = 4
SECTION_ORBITAL_ELEMENTS = 5
SECTION_COMPRESSED_CHEBYSHEV = 6  # LEB2: per-body compressed blobs

# LEB2 format constants
LEB2_MAGIC = b"LEB2"
LEB2_VERSION = 2  # v2: chunked per-body compression (10-year chunks)
LEB2_VERSION_V1 = 1  # v1: monolithic per-body compression (legacy)
COMPRESSION_NONE = 0
COMPRESSION_ZSTD_TRUNC_SHUFFLE = (
    1  # mantissa truncation + coeff-major reorder + byte shuffle + zstd
)

# LEB2 v2 chunk index constants (10-year temporal chunks per body)
CHUNK_INTERVAL_DAYS = 10 * 365.25  # 10 years in days
CHUNK_INDEX_HEADER_FMT = "<IIfI"  # chunk_count, reserved, chunk_interval_days, reserved2
CHUNK_INDEX_HEADER_SIZE = struct.calcsize(CHUNK_INDEX_HEADER_FMT)  # 16
CHUNK_ENTRY_FMT = "<ddQQQII"  # jd_start, jd_end, blob_offset, compressed_size, uncompressed_size, segment_start, segment_count
CHUNK_ENTRY_SIZE = struct.calcsize(CHUNK_ENTRY_FMT)  # 48

# Struct formats (little-endian, no padding with '<' prefix)
HEADER_FMT = "<4sIIIdddI20s"
HEADER_SIZE = struct.calcsize(HEADER_FMT)  # 64

SECTION_DIR_FMT = "<IIQQ"
SECTION_DIR_SIZE = struct.calcsize(SECTION_DIR_FMT)  # 24

BODY_ENTRY_FMT = "<iIIdddIIQ"
BODY_ENTRY_SIZE = struct.calcsize(BODY_ENTRY_FMT)  # 52

# LEB2 compressed body entry: LEB1 fields + compressed_size + uncompressed_size
COMPRESSED_BODY_ENTRY_FMT = "<iIIdddIIQQQ"
COMPRESSED_BODY_ENTRY_SIZE = struct.calcsize(COMPRESSED_BODY_ENTRY_FMT)  # 68

STAR_ENTRY_FMT = "<iddddddd4s"
STAR_ENTRY_SIZE = struct.calcsize(STAR_ENTRY_FMT)  # 64

# Nutation section header format
NUTATION_HEADER_FMT = "<dddIIII"
NUTATION_HEADER_SIZE = struct.calcsize(NUTATION_HEADER_FMT)  # 40

# Delta-T entry format (jd, delta_t_days)
DELTA_T_ENTRY_FMT = "<dd"
DELTA_T_ENTRY_SIZE = struct.calcsize(DELTA_T_ENTRY_FMT)  # 16

# Delta-T section header
DELTA_T_HEADER_FMT = "<II"
DELTA_T_HEADER_SIZE = struct.calcsize(DELTA_T_HEADER_FMT)  # 8


# =============================================================================
# DATACLASSES
# =============================================================================


@dataclass
class FileHeader:
    """LEB file header (64 bytes)."""

    magic: bytes
    version: int
    section_count: int
    body_count: int
    jd_start: float
    jd_end: float
    generation_epoch: float
    flags: int


@dataclass
class SectionEntry:
    """Section directory entry (24 bytes)."""

    section_id: int
    offset: int  # Byte offset from file start
    size: int  # Section size in bytes


@dataclass
class BodyEntry:
    """Body index entry (52 bytes)."""

    body_id: int
    coord_type: int  # COORD_ICRS_BARY, COORD_ECLIPTIC, COORD_HELIO_ECL
    segment_count: int
    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int  # 3 for xyz/lonlatdist
    data_offset: int  # Absolute byte offset to first coefficient


@dataclass
class StarEntry:
    """Fixed star catalog entry (64 bytes)."""

    star_id: int
    ra_j2000: float  # degrees
    dec_j2000: float  # degrees
    pm_ra: float  # deg/yr (includes cos(dec) factor)
    pm_dec: float  # deg/yr
    parallax: float  # arcsec
    rv: float  # km/s
    magnitude: float


@dataclass
class CompressedBodyEntry:
    """LEB2 compressed body index entry (68 bytes).

    Same fields as BodyEntry plus compressed/uncompressed sizes.
    """

    body_id: int
    coord_type: int
    segment_count: int
    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int
    data_offset: int  # Offset to compressed blob
    compressed_size: int
    uncompressed_size: int


@dataclass
class ChunkEntry:
    """LEB2 v2 chunk index entry (48 bytes).

    Each chunk covers a temporal slice (typically 10 years) of a body's
    Chebyshev coefficients.  Only the chunk containing the requested JD
    is decompressed, avoiding the full-body decompression penalty of v1.
    """

    jd_start: float  # Start JD of this chunk
    jd_end: float  # End JD of this chunk
    blob_offset: int  # Absolute byte offset to compressed blob
    compressed_size: int  # Bytes in compressed blob
    uncompressed_size: int  # Original bytes before compression
    segment_start: int  # First segment index in this chunk
    segment_count: int  # Number of segments in this chunk


@dataclass
class NutationHeader:
    """Nutation section header (40 bytes)."""

    jd_start: float
    jd_end: float
    interval_days: float
    degree: int
    components: int  # 2 (dpsi, deps)
    segment_count: int
    reserved: int  # Padding


# =============================================================================
# BODY PARAMETER TABLE
# =============================================================================
# Maps body_id -> (interval_days, degree, coord_type, components)
# This is the single source of truth for Chebyshev parameters.

BODY_PARAMS: dict[int, tuple[float, int, int, int]] = {
    # SUN through MARS, EARTH: ICRS barycentric planet center
    # Parameters: (interval_days, degree, coord_type, components)
    # Smooth barycentric trajectories fit Chebyshev polynomials well.
    # Runtime pipeline applies gravitational deflection + SR aberration
    # for <0.001" total error vs Skyfield apparent().
    0: (32, 13, COORD_ICRS_BARY, 3),  # SUN       — <0.001"
    1: (
        4,
        13,
        COORD_ICRS_BARY,
        3,
    ),  # MOON      — <0.001" (margin 2.8x, DO NOT CHANGE)
    2: (16, 15, COORD_ICRS_BARY, 3),  # MERCURY   — <0.001"
    3: (16, 13, COORD_ICRS_BARY, 3),  # VENUS     — <0.001"
    4: (16, 13, COORD_ICRS_BARY, 3),  # MARS      — <0.001"
    # Outer planets: ICRS system barycenter — COB correction applied at runtime.
    5: (32, 13, COORD_ICRS_BARY_SYSTEM, 3),  # JUPITER — <0.001"
    6: (32, 13, COORD_ICRS_BARY_SYSTEM, 3),  # SATURN  — <0.001"
    7: (64, 13, COORD_ICRS_BARY_SYSTEM, 3),  # URANUS  — <0.001"
    8: (64, 13, COORD_ICRS_BARY_SYSTEM, 3),  # NEPTUNE — <0.001"
    9: (64, 11, COORD_ICRS_BARY_SYSTEM, 3),  # PLUTO   — 0.0005" fit error
    14: (4, 13, COORD_ICRS_BARY, 3),  # EARTH     — <0.001" (mirror of Moon)
    # Lunar nodes/Lilith: ecliptic direct
    # These have complex oscillations from lunar perturbations/precession —
    # keep conservative params. Will optimize after dedicated sweep.
    10: (8, 13, COORD_ECLIPTIC, 3),  # MEAN_NODE  — keep current
    11: (8, 13, COORD_ECLIPTIC, 3),  # TRUE_NODE  — keep current
    12: (8, 13, COORD_ECLIPTIC, 3),  # MEAN_APOG  — keep current
    13: (4, 15, COORD_ECLIPTIC, 3),  # OSCU_APOG  — keep current
    21: (
        1,
        17,
        COORD_ECLIPTIC,
        3,
    ),  # INTP_APOG  — 1d interval for rapid oscillations
    22: (
        1,
        17,
        COORD_ECLIPTIC,
        3,
    ),  # INTP_PERG  — 1d interval for rapid oscillations
    # Main asteroids: ICRS barycentric — keep current (16d fails: 15-321")
    15: (8, 13, COORD_ICRS_BARY, 3),  # CHIRON
    17: (8, 13, COORD_ICRS_BARY, 3),  # CERES
    18: (8, 13, COORD_ICRS_BARY, 3),  # PALLAS
    19: (8, 13, COORD_ICRS_BARY, 3),  # JUNO
    20: (8, 13, COORD_ICRS_BARY, 3),  # VESTA
    # Uranian hypotheticals: heliocentric ecliptic — pure analytical, ~0" error
    40: (256, 7, COORD_HELIO_ECL, 3),  # CUPIDO
    41: (256, 7, COORD_HELIO_ECL, 3),  # HADES
    42: (256, 7, COORD_HELIO_ECL, 3),  # ZEUS
    43: (256, 7, COORD_HELIO_ECL, 3),  # KRONOS
    44: (256, 7, COORD_HELIO_ECL, 3),  # APOLLON
    45: (256, 7, COORD_HELIO_ECL, 3),  # ADMETOS
    46: (256, 7, COORD_HELIO_ECL, 3),  # VULKANUS
    47: (256, 7, COORD_HELIO_ECL, 3),  # POSEIDON
    48: (256, 7, COORD_HELIO_ECL, 3),  # ISIS
}


# =============================================================================
# SERIALIZATION HELPERS
# =============================================================================
# All use struct.pack_into / struct.unpack_from for zero-copy mmap compatibility.


def write_header(buf: bytearray, header: FileHeader) -> None:
    """Pack a FileHeader into a bytearray at offset 0."""
    struct.pack_into(
        HEADER_FMT,
        buf,
        0,
        header.magic,
        header.version,
        header.section_count,
        header.body_count,
        header.jd_start,
        header.jd_end,
        header.generation_epoch,
        header.flags,
        b"\x00" * 20,
    )


def read_header(data: Any, offset: int = 0) -> FileHeader:
    """Unpack a FileHeader from bytes."""
    fields = struct.unpack_from(HEADER_FMT, data, offset)
    return FileHeader(
        magic=fields[0],
        version=fields[1],
        section_count=fields[2],
        body_count=fields[3],
        jd_start=fields[4],
        jd_end=fields[5],
        generation_epoch=fields[6],
        flags=fields[7],
    )


def write_section_dir(buf: bytearray, offset: int, entry: SectionEntry) -> None:
    """Pack a SectionEntry into a bytearray at the given offset."""
    struct.pack_into(
        SECTION_DIR_FMT,
        buf,
        offset,
        entry.section_id,
        0,  # reserved padding
        entry.offset,
        entry.size,
    )


def read_section_dir(data: Any, offset: int) -> SectionEntry:
    """Unpack a SectionEntry from bytes at the given offset."""
    fields = struct.unpack_from(SECTION_DIR_FMT, data, offset)
    return SectionEntry(
        section_id=fields[0],
        # fields[1] is reserved padding
        offset=fields[2],
        size=fields[3],
    )


def write_body_entry(buf: bytearray, offset: int, entry: BodyEntry) -> None:
    """Pack a BodyEntry into a bytearray at the given offset."""
    struct.pack_into(
        BODY_ENTRY_FMT,
        buf,
        offset,
        entry.body_id,
        entry.coord_type,
        entry.segment_count,
        entry.jd_start,
        entry.jd_end,
        entry.interval_days,
        entry.degree,
        entry.components,
        entry.data_offset,
    )


def read_body_entry(data: Any, offset: int) -> BodyEntry:
    """Unpack a BodyEntry from bytes at the given offset."""
    fields = struct.unpack_from(BODY_ENTRY_FMT, data, offset)
    return BodyEntry(
        body_id=fields[0],
        coord_type=fields[1],
        segment_count=fields[2],
        jd_start=fields[3],
        jd_end=fields[4],
        interval_days=fields[5],
        degree=fields[6],
        components=fields[7],
        data_offset=fields[8],
    )


def write_star_entry(buf: bytearray, offset: int, entry: StarEntry) -> None:
    """Pack a StarEntry into a bytearray at the given offset."""
    struct.pack_into(
        STAR_ENTRY_FMT,
        buf,
        offset,
        entry.star_id,
        entry.ra_j2000,
        entry.dec_j2000,
        entry.pm_ra,
        entry.pm_dec,
        entry.parallax,
        entry.rv,
        entry.magnitude,
        b"\x00" * 4,
    )


def read_star_entry(data: Any, offset: int) -> StarEntry:
    """Unpack a StarEntry from bytes at the given offset."""
    fields = struct.unpack_from(STAR_ENTRY_FMT, data, offset)
    return StarEntry(
        star_id=fields[0],
        ra_j2000=fields[1],
        dec_j2000=fields[2],
        pm_ra=fields[3],
        pm_dec=fields[4],
        parallax=fields[5],
        rv=fields[6],
        magnitude=fields[7],
    )


def write_nutation_header(buf: bytearray, offset: int, header: NutationHeader) -> None:
    """Pack a NutationHeader into a bytearray at the given offset."""
    struct.pack_into(
        NUTATION_HEADER_FMT,
        buf,
        offset,
        header.jd_start,
        header.jd_end,
        header.interval_days,
        header.degree,
        header.components,
        header.segment_count,
        header.reserved,
    )


def read_nutation_header(data: Any, offset: int) -> NutationHeader:
    """Unpack a NutationHeader from bytes at the given offset."""
    fields = struct.unpack_from(NUTATION_HEADER_FMT, data, offset)
    return NutationHeader(
        jd_start=fields[0],
        jd_end=fields[1],
        interval_days=fields[2],
        degree=fields[3],
        components=fields[4],
        segment_count=fields[5],
        reserved=fields[6],
    )


def write_compressed_body_entry(
    buf: bytearray, offset: int, entry: CompressedBodyEntry
) -> None:
    """Pack a CompressedBodyEntry into a bytearray at the given offset."""
    struct.pack_into(
        COMPRESSED_BODY_ENTRY_FMT,
        buf,
        offset,
        entry.body_id,
        entry.coord_type,
        entry.segment_count,
        entry.jd_start,
        entry.jd_end,
        entry.interval_days,
        entry.degree,
        entry.components,
        entry.data_offset,
        entry.compressed_size,
        entry.uncompressed_size,
    )


def read_compressed_body_entry(data: Any, offset: int) -> CompressedBodyEntry:
    """Unpack a CompressedBodyEntry from bytes at the given offset."""
    fields = struct.unpack_from(COMPRESSED_BODY_ENTRY_FMT, data, offset)
    return CompressedBodyEntry(
        body_id=fields[0],
        coord_type=fields[1],
        segment_count=fields[2],
        jd_start=fields[3],
        jd_end=fields[4],
        interval_days=fields[5],
        degree=fields[6],
        components=fields[7],
        data_offset=fields[8],
        compressed_size=fields[9],
        uncompressed_size=fields[10],
    )


def write_chunk_index_header(
    buf: bytearray, offset: int, chunk_count: int, chunk_interval_days: float
) -> None:
    """Pack a chunk index header into a bytearray."""
    struct.pack_into(
        CHUNK_INDEX_HEADER_FMT,
        buf,
        offset,
        chunk_count,
        0,  # reserved
        chunk_interval_days,
        0,  # reserved2
    )


def read_chunk_index_header(
    data: Any, offset: int
) -> tuple:
    """Unpack a chunk index header. Returns (chunk_count, chunk_interval_days)."""
    fields = struct.unpack_from(CHUNK_INDEX_HEADER_FMT, data, offset)
    return fields[0], fields[2]  # chunk_count, chunk_interval_days


def write_chunk_entry(buf: bytearray, offset: int, entry: ChunkEntry) -> None:
    """Pack a ChunkEntry into a bytearray."""
    struct.pack_into(
        CHUNK_ENTRY_FMT,
        buf,
        offset,
        entry.jd_start,
        entry.jd_end,
        entry.blob_offset,
        entry.compressed_size,
        entry.uncompressed_size,
        entry.segment_start,
        entry.segment_count,
    )


def read_chunk_entry(data: Any, offset: int) -> ChunkEntry:
    """Unpack a ChunkEntry from bytes."""
    fields = struct.unpack_from(CHUNK_ENTRY_FMT, data, offset)
    return ChunkEntry(
        jd_start=fields[0],
        jd_end=fields[1],
        blob_offset=fields[2],
        compressed_size=fields[3],
        uncompressed_size=fields[4],
        segment_start=fields[5],
        segment_count=fields[6],
    )


def segment_byte_size(degree: int, components: int) -> int:
    """Calculate the byte size of a single Chebyshev segment.

    Args:
        degree: Polynomial degree.
        components: Number of components (3 for xyz or lonlatdist).

    Returns:
        Size in bytes of one segment's coefficient data.
    """
    return components * (degree + 1) * 8  # 8 bytes per f64


# ---------------------------------------------------------------------------
# mmap warm helper
# ---------------------------------------------------------------------------

try:
    _PAGE_SIZE = os.sysconf("SC_PAGE_SIZE")
except (AttributeError, OSError, ValueError):  # pragma: no cover - POSIX
    # sysconf is always present on the supported platforms (macOS/Linux);
    # the 4096 fallback is a defensive default for exotic environments.
    _PAGE_SIZE = 4096
_PAGE_MASK = ~(_PAGE_SIZE - 1)


def _madvise_ranges(mm: Any, ranges: list[tuple[int, int]]) -> None:
    """Call madvise(MADV_WILLNEED) on merged, page-aligned byte ranges.

    Args:
        mm: An open mmap object.
        ranges: List of (offset, size) byte ranges to pre-fault.
    """
    import mmap as _mmap

    if not ranges:
        return

    ranges.sort()
    merged: list[tuple[int, int]] = [ranges[0]]
    for offset, size in ranges[1:]:
        prev_off, prev_size = merged[-1]
        if offset <= prev_off + prev_size + _PAGE_SIZE:
            merged[-1] = (prev_off, max(prev_size, offset + size - prev_off))
        else:
            merged.append((offset, size))

    mm_len = len(mm)
    for offset, size in merged:
        start_page = offset & _PAGE_MASK
        end_page = (offset + size + _PAGE_SIZE - 1) & _PAGE_MASK
        length = min(end_page, mm_len) - start_page
        if length <= 0:
            continue
        try:
            mm.madvise(_mmap.MADV_WILLNEED, start_page, length)
        except (OSError, AttributeError, ValueError):
            pass


def _madvise_dontneed(mm: Any) -> None:
    """Advise the kernel that all mmap pages can be reclaimed.

    Args:
        mm: An open mmap object (or ``None`` for a no-op).
    """
    import mmap as _mmap

    if mm is None:
        return
    try:
        mm.madvise(_mmap.MADV_DONTNEED)
    except (OSError, AttributeError, ValueError):
        pass
