#!/usr/bin/env python3
"""
LEB2 compressed ephemeris file generator.

Provenance:
    Project converter/generator over registered LEB/JPL/IAU inputs. Mantissa
    truncation follows IEEE-754 binary64; coefficient ordering and byte shuffle
    are documented reversible transforms; final compression uses Zstandard.
    Per-body bit counts are derived from explicit error targets and verified
    against the uncompressed values. No hidden reconstruction table is used.

Two modes of operation:

1. Convert: Read an existing LEB1 file and produce a compressed LEB2 file.
   python scripts/generate_leb2.py convert data/leb/ephemeris_base.leb -o core_base.leb --group core --tier base

2. Generate: Generate from scratch via Skyfield + Chebyshev fitting, then compress.
   python scripts/generate_leb2.py generate --tier base --group core -o core_base.leb

The LEB2 format uses error-bounded lossy compression (mantissa truncation +
coefficient-major reorder + byte shuffle + zstd) achieving 5-15x compression
per body. Coefficient targets use each body's native stored units: AU for ICRS
Cartesian data and degrees/degrees/AU for ecliptic data.

Body groups:
  core       : Sun, Moon, Mercury-Pluto, Earth, Mean/True Node, Mean Apogee (14 bodies)
  asteroids  : Chiron, Ceres, Pallas, Juno, Vesta (5 bodies)
  exotics    : Centaurs, TNOs, and NEAs from the canonical exotic registry
  apogee     : Osculating Apogee, Interpolated Apogee/Perigee (3 bodies)
"""

from __future__ import annotations

import argparse
import mmap
import os
import sys
import time
from typing import Iterable, Optional

import numpy as np

# Allow importing from the library
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from libephemeris.leb_compression import (
    BODY_TARGET_AU,
    DEFAULT_TARGET_AU,
    compress_body_chunked,
    compute_mantissa_bits,
)
from libephemeris.leb_format import (
    BODY_ENTRY_SIZE,
    CHUNK_ENTRY_SIZE,
    CHUNK_INDEX_HEADER_SIZE,
    CHUNK_INTERVAL_DAYS,
    COMPRESSED_BODY_ENTRY_SIZE,
    COMPRESSION_ZSTD_TRUNC_SHUFFLE,
    COORD_ECLIPTIC,
    COORD_GEO_ECLIPTIC,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
    HEADER_SIZE,
    LEB2_MAGIC,
    LEB2_VERSION,
    LEB2_VERSION_V1,
    MAGIC,
    SECTION_BODY_INDEX,
    SECTION_COMPRESSED_CHEBYSHEV,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    SECTION_STARS,
    ChunkEntry,
    CompressedBodyEntry,
    FileHeader,
    SectionEntry,
    read_body_entry,
    read_header,
    read_section_dir,
    segment_byte_size,
    write_chunk_entry,
    write_chunk_index_header,
    write_compressed_body_entry,
    write_header,
    write_section_dir,
)
from libephemeris.leb_groups import LEB2_GROUPS as CANONICAL_LEB2_GROUPS
from libephemeris.exotic_bodies import (
    EXOTIC_EXTENDED_IDS,
    EXOTIC_IDS,
    name_map as _exotic_names,
)

# =============================================================================
# BODY GROUPS
# =============================================================================

LEB2_GROUPS = {
    "core": sorted([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 10, 11, 12]),
    "asteroids": sorted([15, 17, 18, 19, 20]),
    "exotics": list(EXOTIC_IDS),  # centaurs/TNOs/NEAs — registry source of truth
    "apogee": sorted([13, 21, 22]),
}
assert tuple(LEB2_GROUPS) == CANONICAL_LEB2_GROUPS

_ALLOWED_GROUP_INVENTORIES: dict[str, tuple[frozenset[int], ...]] = {
    "core": (frozenset(LEB2_GROUPS["core"]),),
    "asteroids": (frozenset(LEB2_GROUPS["asteroids"]),),
    "exotics": (frozenset(EXOTIC_IDS),),
    "apogee": (frozenset(LEB2_GROUPS["apogee"]),),
}
_ANGULAR_COORD_TYPES = {COORD_ECLIPTIC, COORD_GEO_ECLIPTIC}
_CARTESIAN_COORD_TYPES = {COORD_ICRS_BARY, COORD_ICRS_BARY_SYSTEM}


def _group_inventory_error(
    group: str,
    body_ids: Iterable[int],
    expected_tier: Optional[str] = None,
) -> Optional[str]:
    """Return an error for a non-canonical named-group body inventory.

    The generic and base/medium exotic inventories contain all 31 registered
    bodies. Only the explicitly declared extended tier may use the 23-body
    non-NEA subset.
    """
    if expected_tier not in (None, "base", "medium", "extended"):
        raise ValueError(f"unknown LEB tier {expected_tier!r}")
    actual = frozenset(body_ids)
    if group == "exotics" and expected_tier == "extended":
        allowed = (frozenset(EXOTIC_EXTENDED_IDS),)
    else:
        allowed = _ALLOWED_GROUP_INVENTORIES[group]
    if actual in allowed:
        return None
    expected = " or ".join(str(sorted(inventory)) for inventory in allowed)
    return (
        f"invalid {group!r} body inventory: got {sorted(actual)}; "
        f"expected exactly {expected}"
    )


BODY_NAMES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
    10: "Mean Node",
    11: "True Node",
    12: "Mean Apogee",
    13: "Oscu Apogee",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "Interp Apogee",
    22: "Interp Perigee",
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Transpluto",
}
BODY_NAMES.update(_exotic_names())  # exotic minor bodies — registry labels


# =============================================================================
# LEB1 READER (minimal, for conversion)
# =============================================================================


class LEB1Source:
    """Reads an LEB1 file and extracts raw data for conversion to LEB2."""

    def __init__(self, path: str):
        self.path = path
        self._file = open(path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._parse()

    def _parse(self):
        self.header = read_header(self._mm, 0)
        if self.header.magic != MAGIC:
            raise ValueError(f"Not a LEB1 file: magic={self.header.magic!r}")

        self.sections = {}
        for i in range(self.header.section_count):
            sec = read_section_dir(self._mm, HEADER_SIZE + i * SECTION_DIR_SIZE)
            self.sections[sec.section_id] = sec

        self.bodies = {}
        if SECTION_BODY_INDEX in self.sections:
            sec = self.sections[SECTION_BODY_INDEX]
            for i in range(self.header.body_count):
                entry = read_body_entry(self._mm, sec.offset + i * BODY_ENTRY_SIZE)
                self.bodies[entry.body_id] = entry

    def get_body_coefficients(self, body_id: int) -> bytes:
        """Extract raw coefficient bytes for a body."""
        entry = self.bodies[body_id]
        seg_size = segment_byte_size(entry.degree, entry.components)
        raw_size = seg_size * entry.segment_count
        return bytes(self._mm[entry.data_offset : entry.data_offset + raw_size])

    def get_section_data(self, section_id: int) -> Optional[bytes]:
        """Extract raw bytes for a section (nutation, delta-t, stars)."""
        if section_id not in self.sections:
            return None
        sec = self.sections[section_id]
        return bytes(self._mm[sec.offset : sec.offset + sec.size])

    def close(self):
        if self._mm:
            self._mm.close()
        if self._file:
            self._file.close()


# =============================================================================
# LEB2 WRITER
# =============================================================================


def write_leb2(
    output: str,
    body_entries: list,  # list of (BodyEntry, compressed_bytes, uncompressed_size)
    nutation_data: Optional[bytes],
    delta_t_data: Optional[bytes],
    star_data: Optional[bytes],
    jd_start: float,
    jd_end: float,
    generation_epoch: float,
    verbose: bool = True,
) -> None:
    """Write a LEB2 v1 (monolithic) file from compressed body data.

    Each body is stored as a single compressed blob (no per-body chunk
    index).  The header is stamped with ``LEB2_VERSION_V1`` so ``open_leb``
    routes reads through the v1 monolithic decode path.  ``write_leb2_chunked``
    produces the v2 chunked format used by the current conversion pipeline.
    """

    n_bodies = len(body_entries)

    # Count active sections
    section_list = [SECTION_BODY_INDEX, SECTION_COMPRESSED_CHEBYSHEV]
    if nutation_data:
        section_list.append(SECTION_NUTATION)
    if delta_t_data:
        section_list.append(SECTION_DELTA_T)
    if star_data:
        section_list.append(SECTION_STARS)
    n_sections = len(section_list)

    # Calculate layout
    body_index_offset = HEADER_SIZE + n_sections * SECTION_DIR_SIZE
    body_index_size = n_bodies * COMPRESSED_BODY_ENTRY_SIZE

    cheb_offset = body_index_offset + body_index_size
    total_cheb = sum(len(blob) for _, blob, _ in body_entries)

    current = cheb_offset + total_cheb

    nut_offset = current
    if nutation_data:
        current += len(nutation_data)

    dt_offset = current
    if delta_t_data:
        current += len(delta_t_data)

    star_offset = current
    if star_data:
        current += len(star_data)

    total_size = current

    # Allocate buffer
    buf = bytearray(total_size)

    # Write header. This monolithic layout stores one compressed blob per
    # body (no chunk index), which is the LEB2 v1 format — the reader keys
    # off the version to pick the v1 (monolithic) vs v2 (chunked) decode
    # path, so stamping v2 here would make the file unreadable.
    file_hdr = FileHeader(
        magic=LEB2_MAGIC,
        version=LEB2_VERSION_V1,
        section_count=n_sections,
        body_count=n_bodies,
        jd_start=jd_start,
        jd_end=jd_end,
        generation_epoch=generation_epoch,
        flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
    )
    write_header(buf, file_hdr)

    # Write section directory
    sec_entries = [
        SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
        SectionEntry(SECTION_COMPRESSED_CHEBYSHEV, cheb_offset, total_cheb),
    ]
    if nutation_data:
        sec_entries.append(
            SectionEntry(SECTION_NUTATION, nut_offset, len(nutation_data))
        )
    if delta_t_data:
        sec_entries.append(SectionEntry(SECTION_DELTA_T, dt_offset, len(delta_t_data)))
    if star_data:
        sec_entries.append(SectionEntry(SECTION_STARS, star_offset, len(star_data)))

    for i, se in enumerate(sec_entries):
        write_section_dir(buf, HEADER_SIZE + i * SECTION_DIR_SIZE, se)

    # Write body index and compressed blobs
    blob_offset = cheb_offset
    for idx, (body, blob, uncompressed_size) in enumerate(body_entries):
        cbe = CompressedBodyEntry(
            body_id=body.body_id,
            coord_type=body.coord_type,
            segment_count=body.segment_count,
            jd_start=body.jd_start,
            jd_end=body.jd_end,
            interval_days=body.interval_days,
            degree=body.degree,
            components=body.components,
            data_offset=blob_offset,
            compressed_size=len(blob),
            uncompressed_size=uncompressed_size,
        )
        write_compressed_body_entry(
            buf, body_index_offset + idx * COMPRESSED_BODY_ENTRY_SIZE, cbe
        )
        buf[blob_offset : blob_offset + len(blob)] = blob
        blob_offset += len(blob)

    # Write auxiliary sections
    if nutation_data:
        buf[nut_offset : nut_offset + len(nutation_data)] = nutation_data
    if delta_t_data:
        buf[dt_offset : dt_offset + len(delta_t_data)] = delta_t_data
    if star_data:
        buf[star_offset : star_offset + len(star_data)] = star_data

    # Write to disk
    os.makedirs(os.path.dirname(os.path.abspath(output)) or ".", exist_ok=True)
    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print(f"\n  Output: {output}")
        print(f"  File size: {total_size / 1e6:.2f} MB ({total_size:,} bytes)")


# =============================================================================
# LEB2 v2 CHUNKED WRITER
# =============================================================================


def write_leb2_chunked(
    output: str,
    body_chunks: list,
    # list of (BodyEntry, list[(compressed_bytes, uncompressed_size)], total_uncompressed)
    nutation_data: Optional[bytes],
    delta_t_data: Optional[bytes],
    star_data: Optional[bytes],
    jd_start: float,
    jd_end: float,
    generation_epoch: float,
    chunk_interval_days: float = CHUNK_INTERVAL_DAYS,
    verbose: bool = True,
) -> None:
    """Write a LEB2 v2 file with chunked per-body compression.

    Each body's data is split into temporal chunks (default 10 years).
    Each chunk is compressed independently, enabling per-chunk decompression
    at read time for near-instant cold-start performance.
    """
    n_bodies = len(body_chunks)

    # Count sections
    section_list = [SECTION_BODY_INDEX, SECTION_COMPRESSED_CHEBYSHEV]
    if nutation_data:
        section_list.append(SECTION_NUTATION)
    if delta_t_data:
        section_list.append(SECTION_DELTA_T)
    if star_data:
        section_list.append(SECTION_STARS)
    n_sections = len(section_list)

    # Calculate layout: header + section dir + body index
    body_index_offset = HEADER_SIZE + n_sections * SECTION_DIR_SIZE
    body_index_size = n_bodies * COMPRESSED_BODY_ENTRY_SIZE

    # Chebyshev data: for each body = chunk_index_header + chunk_entries + chunk_blobs
    cheb_offset = body_index_offset + body_index_size
    cheb_total = 0
    for body, chunks, _ in body_chunks:
        n_chunks = len(chunks)
        chunk_index_size = CHUNK_INDEX_HEADER_SIZE + n_chunks * CHUNK_ENTRY_SIZE
        blobs_size = sum(len(blob) for blob, _ in chunks)
        cheb_total += chunk_index_size + blobs_size

    current = cheb_offset + cheb_total

    nut_offset = current
    if nutation_data:
        current += len(nutation_data)
    dt_offset = current
    if delta_t_data:
        current += len(delta_t_data)
    star_offset = current
    if star_data:
        current += len(star_data)

    total_size = current
    buf = bytearray(total_size)

    # Write header (version 2)
    file_hdr = FileHeader(
        magic=LEB2_MAGIC,
        version=LEB2_VERSION,
        section_count=n_sections,
        body_count=n_bodies,
        jd_start=jd_start,
        jd_end=jd_end,
        generation_epoch=generation_epoch,
        flags=COMPRESSION_ZSTD_TRUNC_SHUFFLE,
    )
    write_header(buf, file_hdr)

    # Write section directory
    sec_entries = [
        SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
        SectionEntry(SECTION_COMPRESSED_CHEBYSHEV, cheb_offset, cheb_total),
    ]
    if nutation_data:
        sec_entries.append(
            SectionEntry(SECTION_NUTATION, nut_offset, len(nutation_data))
        )
    if delta_t_data:
        sec_entries.append(SectionEntry(SECTION_DELTA_T, dt_offset, len(delta_t_data)))
    if star_data:
        sec_entries.append(SectionEntry(SECTION_STARS, star_offset, len(star_data)))
    for i, se in enumerate(sec_entries):
        write_section_dir(buf, HEADER_SIZE + i * SECTION_DIR_SIZE, se)

    # Write body index and chunked data
    data_cursor = cheb_offset
    for idx, (body, chunks, total_uncomp) in enumerate(body_chunks):
        n_chunks = len(chunks)
        chunk_index_size = CHUNK_INDEX_HEADER_SIZE + n_chunks * CHUNK_ENTRY_SIZE
        blobs_size = sum(len(blob) for blob, _ in chunks)
        body_data_size = chunk_index_size + blobs_size

        # Body entry: data_offset points to the chunk index header
        cbe = CompressedBodyEntry(
            body_id=body.body_id,
            coord_type=body.coord_type,
            segment_count=body.segment_count,
            jd_start=body.jd_start,
            jd_end=body.jd_end,
            interval_days=body.interval_days,
            degree=body.degree,
            components=body.components,
            data_offset=data_cursor,
            compressed_size=body_data_size,
            uncompressed_size=total_uncomp,
        )
        write_compressed_body_entry(
            buf, body_index_offset + idx * COMPRESSED_BODY_ENTRY_SIZE, cbe
        )

        # Write chunk index header
        write_chunk_index_header(buf, data_cursor, n_chunks, chunk_interval_days)
        ci_cursor = data_cursor + CHUNK_INDEX_HEADER_SIZE

        # Blob data starts after the chunk index
        blob_base = data_cursor + chunk_index_size
        blob_cursor = blob_base

        # Build chunk entries
        seg_bytes = segment_byte_size(body.degree, body.components)
        seg_cursor = 0  # running segment start

        for ci, (blob, uncomp_size) in enumerate(chunks):
            chunk_seg_count = uncomp_size // seg_bytes
            chunk_jd_start = body.jd_start + seg_cursor * body.interval_days
            chunk_jd_end = chunk_jd_start + chunk_seg_count * body.interval_days

            chunk_entry = ChunkEntry(
                jd_start=chunk_jd_start,
                jd_end=chunk_jd_end,
                blob_offset=blob_cursor,  # absolute offset
                compressed_size=len(blob),
                uncompressed_size=uncomp_size,
                segment_start=seg_cursor,
                segment_count=chunk_seg_count,
            )
            write_chunk_entry(buf, ci_cursor + ci * CHUNK_ENTRY_SIZE, chunk_entry)

            # Write blob
            buf[blob_cursor : blob_cursor + len(blob)] = blob
            blob_cursor += len(blob)
            seg_cursor += chunk_seg_count

        data_cursor += body_data_size

    # Write auxiliary sections
    if nutation_data:
        buf[nut_offset : nut_offset + len(nutation_data)] = nutation_data
    if delta_t_data:
        buf[dt_offset : dt_offset + len(delta_t_data)] = delta_t_data
    if star_data:
        buf[star_offset : star_offset + len(star_data)] = star_data

    # Flush
    os.makedirs(os.path.dirname(os.path.abspath(output)) or ".", exist_ok=True)
    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print(f"\n  Output: {output}")
        print(f"  File size: {total_size / 1e6:.2f} MB ({total_size:,} bytes)")
        print(
            f"  Format: LEB2 v2 (chunked, {chunk_interval_days / 365.25:.0f}-year chunks)"
        )


# =============================================================================
# CONVERT MODE: LEB1 -> LEB2 (v2 chunked)
# =============================================================================


def convert_leb1_to_leb2(
    input_path: str,
    output_path: str,
    group: Optional[str] = None,
    expected_tier: Optional[str] = None,
    target_precision: float = DEFAULT_TARGET_AU,
    verbose: bool = True,
) -> None:
    """Convert an existing LEB1 file to LEB2 v2 format (chunked).

    Args:
        input_path: Path to source LEB1 file.
        output_path: Path to output LEB2 file.
        group: Optional body group filter (core/asteroids/exotics/apogee).
            Named companion groups contain body channels only; shared
            nutation, Delta-T, and star-catalog sections belong to ``core``.
        expected_tier: Optional tier used to authenticate tier-specific inventory.
        target_precision: Numeric coefficient target in native stored units.
        verbose: Print progress.
    """
    if verbose:
        chunk_years = CHUNK_INTERVAL_DAYS / 365.25
        print(
            f"Converting {input_path} -> {output_path} (v2, {chunk_years:.0f}-year chunks)"
        )
        if group:
            print(f"  Group: {group}")

    t0 = time.time()
    src = LEB1Source(input_path)

    # Select bodies
    if group:
        body_ids = [bid for bid in LEB2_GROUPS[group] if bid in src.bodies]
    else:
        body_ids = sorted(src.bodies.keys())
    if group and not body_ids:
        src.close()
        raise ValueError(
            f"source LEB1 contains no bodies for requested group {group!r}"
        )
    if group and (
        inventory_error := _group_inventory_error(group, body_ids, expected_tier)
    ):
        src.close()
        raise ValueError(f"source LEB1 has {inventory_error}")

    if verbose:
        print(f"  Bodies: {len(body_ids)}")
        input_size = os.path.getsize(input_path)
        print(f"  Input size: {input_size / 1e6:.1f} MB")

    # Compress each body into chunks
    body_chunks = []
    total_raw = 0
    total_compressed = 0

    if verbose:
        print(
            f"\n  {'Body':<16s}  {'Raw KB':>8s}  {'Comp KB':>8s}  {'Ratio':>6s}  {'Chunks':>6s}"
        )
        print(f"  {'-' * 52}")

    for bid in body_ids:
        entry = src.bodies[bid]
        raw = src.get_body_coefficients(bid)
        raw_size = len(raw)

        # Compute mantissa bits (use per-body target if available)
        body_target = BODY_TARGET_AU.get(bid, target_precision)
        deg1 = entry.degree + 1
        arr = np.frombuffer(raw, dtype=np.float64).reshape(
            entry.segment_count, entry.components, deg1
        )
        bits = compute_mantissa_bits(arr, body_target)

        # Calculate segments per chunk
        segs_per_chunk = max(1, int(CHUNK_INTERVAL_DAYS / entry.interval_days))

        # Compress into chunks
        chunks = compress_body_chunked(
            raw,
            entry.segment_count,
            entry.degree,
            entry.components,
            bits,
            segs_per_chunk,
        )

        comp_size = sum(len(blob) for blob, _ in chunks)
        body_chunks.append((entry, chunks, raw_size))
        total_raw += raw_size
        total_compressed += comp_size

        if verbose:
            name = BODY_NAMES.get(bid, f"Body {bid}")
            ratio = raw_size / comp_size if comp_size > 0 else 0
            print(
                f"  {name:<16s}  {raw_size / 1024:>8.1f}  {comp_size / 1024:>8.1f}  {ratio:>5.1f}x  {len(chunks):>6d}"
            )

    if verbose:
        print(f"  {'-' * 52}")
        ratio = total_raw / total_compressed if total_compressed > 0 else 0
        print(
            f"  {'TOTAL':<16s}  {total_raw / 1024:>8.1f}  {total_compressed / 1024:>8.1f}  {ratio:>5.1f}x"
        )

    # Shared runtime tables belong to the core artifact.  Copying them into
    # every named companion would duplicate the same large nutation table four
    # times per tier (especially expensive over the full DE441 interval) while
    # CompositeLEBReader already resolves auxiliary data through its primary
    # core reader.  An unfiltered monolithic conversion keeps the historical
    # self-contained behavior.
    include_auxiliary = group in (None, "core")
    nutation_data = (
        src.get_section_data(SECTION_NUTATION) if include_auxiliary else None
    )
    delta_t_data = src.get_section_data(SECTION_DELTA_T) if include_auxiliary else None
    star_data = src.get_section_data(SECTION_STARS) if include_auxiliary else None

    # Write LEB2 v2 (chunked)
    write_leb2_chunked(
        output=output_path,
        body_chunks=body_chunks,
        nutation_data=nutation_data,
        delta_t_data=delta_t_data,
        star_data=star_data,
        jd_start=src.header.jd_start,
        jd_end=src.header.jd_end,
        generation_epoch=src.header.generation_epoch,
        verbose=verbose,
    )

    src.close()

    elapsed = time.time() - t0
    if verbose:
        print(f"  Time: {elapsed:.1f}s")


# =============================================================================
# GENERATE MODE: Skyfield -> Chebyshev -> LEB2
# =============================================================================

TIER_CONFIGS = {
    "base": ("de440s.bsp", 1850, 2150, "ephemeris_base"),
    "medium": ("de440.bsp", 1550, 2650, "ephemeris_medium"),
    "extended": ("de441.bsp", -13200, 17191, "ephemeris_extended"),
}

DEFAULT_LEB_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "leb")


def _year_to_jd(year: int) -> float:
    """Convert year to Julian Day (January 1.0)."""
    if year >= 0:
        return 1721425.5 + 365.25 * year - int(0.01 * year) + int(0.0025 * year)
    a = (14 - 1) // 12
    y = year + 4800 - a
    m = 1 + 12 * a - 3
    return (
        1
        + (153 * m + 2) // 5
        + 365 * y
        + y // 4
        - y // 100
        + y // 400
        - 32045
        + 0.5
        - 0.5
    )


def generate_and_compress(
    output_path: str,
    tier: str,
    group: Optional[str] = None,
    start_year: Optional[int] = None,
    end_year: Optional[int] = None,
    workers: int = 1,
    target_precision: float = DEFAULT_TARGET_AU,
    verbose: bool = True,
) -> None:
    """Generate LEB2 from scratch: Skyfield -> Chebyshev fit -> compress.

    This calls the LEB1 generation pipeline (generate_leb.py's assemble_leb),
    creates a temporary LEB1 file, then converts it to LEB2.

    Args:
        output_path: Path to output LEB2 file.
        tier: Precision tier (base/medium/extended).
        group: Optional body group filter.
        start_year: Override start year.
        end_year: Override end year.
        workers: Parallel workers for Chebyshev fitting.
        target_precision: Numeric coefficient target in native stored units.
        verbose: Print progress.
    """
    import subprocess
    import tempfile

    ephem_file, tier_start, tier_end, tier_name = TIER_CONFIGS[tier]

    # Determine which bodies to generate
    if group:
        body_ids = LEB2_GROUPS[group]
        bodies_arg = ",".join(str(b) for b in body_ids)
    else:
        bodies_arg = None

    # Generate LEB1 to a temporary file
    tmp_leb1 = tempfile.mktemp(suffix=".leb", prefix=f"leb1_{tier_name}_")

    if verbose:
        print(f"Step 1: Generating LEB1 (tier={tier}, group={group or 'all'})...")

    gen_cmd = [
        sys.executable,
        os.path.join(os.path.dirname(__file__), "generate_leb.py"),
        "--tier",
        tier,
        "--output",
        tmp_leb1,
        "--workers",
        str(workers),
    ]
    if start_year is not None:
        gen_cmd += ["--start", str(start_year)]
    if end_year is not None:
        gen_cmd += ["--end", str(end_year)]
    if bodies_arg:
        gen_cmd += ["--bodies", bodies_arg]

    result = subprocess.run(gen_cmd)
    if result.returncode != 0:
        print(
            f"Error: LEB1 generation failed (exit {result.returncode})", file=sys.stderr
        )
        sys.exit(result.returncode)

    if verbose:
        print("\nStep 2: Converting LEB1 -> LEB2...")

    # Convert to LEB2
    convert_leb1_to_leb2(
        input_path=tmp_leb1,
        output_path=output_path,
        group=None,  # already filtered in generation
        target_precision=target_precision,
        verbose=verbose,
    )

    # Cleanup
    if os.path.exists(tmp_leb1):
        os.remove(tmp_leb1)
        if verbose:
            print(f"  Removed temporary LEB1: {os.path.basename(tmp_leb1)}")


# =============================================================================
# VERIFICATION
# =============================================================================


VERIFY_DEFAULT_COEFFICIENT_TARGET = DEFAULT_TARGET_AU


def _verification_tolerance(body_id: int, degree: int) -> float:
    """Return the numeric native-component error bound for a compressed body.

    Compression bounds the absolute error of each Chebyshev coefficient. Since
    ``abs(T_k(x)) <= 1`` on an encoded segment, summing ``degree + 1`` terms
    gives a conservative component bound. The unit follows the stored frame:
    AU for Cartesian ICRS components and degrees/degrees/AU for ecliptic
    longitude, latitude, and distance. Per-body coefficient targets match
    those used by the converter; all other bodies use the named default.
    """
    coefficient_target = BODY_TARGET_AU.get(body_id, VERIFY_DEFAULT_COEFFICIENT_TARGET)
    return coefficient_target * (degree + 1)


def _native_component_errors(
    coord_type: int,
    reference: tuple[float, ...],
    candidate: tuple[float, ...],
) -> tuple[tuple[float, float, float], tuple[str, str, str]]:
    """Return component errors and units without mixing coordinate frames."""
    if coord_type in _ANGULAR_COORD_TYPES:
        longitude = abs((candidate[0] - reference[0] + 180.0) % 360.0 - 180.0)
        return (
            (
                longitude,
                abs(candidate[1] - reference[1]),
                abs(candidate[2] - reference[2]),
            ),
            ("deg", "deg", "AU"),
        )
    if coord_type in _CARTESIAN_COORD_TYPES:
        return (
            (
                abs(candidate[0] - reference[0]),
                abs(candidate[1] - reference[1]),
                abs(candidate[2] - reference[2]),
            ),
            ("AU", "AU", "AU"),
        )
    raise ValueError(f"unsupported LEB coordinate type {coord_type}")


def verify_leb2(
    leb2_path: str,
    reference_leb1: Optional[str] = None,
    n_samples: int = 200,
    verbose: bool = True,
    expected_group: Optional[str] = None,
    expected_tier: Optional[str] = None,
) -> bool:
    """Verify a LEB2 file against a LEB1 reference or run a smoke test.

    Args:
        leb2_path: Path to LEB2 file to verify.
        reference_leb1: Optional LEB1 reference file. If provided, compares
            LEB2 position results against LEB1. If omitted, performs a reader
            smoke test only.
        n_samples: Number of random JDs to sample per body.
        verbose: Print results.
        expected_group: Optional canonical group whose exact inventory is required.
        expected_tier: Optional tier used to authenticate tier-specific inventory.

    Returns:
        True if the inventory and all native stored components pass their bounds.
    """
    from libephemeris.leb2_reader import LEB2Reader

    if n_samples < 1:
        raise ValueError("n_samples must be at least 1")

    reader2 = LEB2Reader(leb2_path)
    all_pass = True

    if not reader2._bodies:
        if verbose:
            print(f"Verifying {leb2_path}: FAIL - no bodies found")
        reader2.close()
        return False
    if expected_group:
        inventory_error = _group_inventory_error(
            expected_group, reader2._bodies, expected_tier
        )
        if inventory_error:
            if verbose:
                print(f"Verifying {leb2_path}: FAIL - {inventory_error}")
            reader2.close()
            return False

    if reference_leb1:
        from libephemeris.leb_reader import LEBReader

        reader1 = LEBReader(reference_leb1)

        leb2_range = reader2.jd_range
        leb1_range = reader1.jd_range
        if not np.all(np.isfinite((*leb2_range, *leb1_range))):
            if verbose:
                print("  File coverage metadata is non-finite: FAIL")
            all_pass = False
        elif leb2_range != leb1_range:
            if verbose:
                print(f"  File coverage {leb2_range} != LEB1 {leb1_range}: FAIL")
            all_pass = False

        if verbose:
            print(f"Verifying {leb2_path} against {reference_leb1}")
            print(f"  Samples per body: {n_samples}")
            print(
                f"\n  {'Body':<16s}  {'Worst component':>18s}  {'Unit':>5s}  "
                f"{'Bound':>10s}  {'Ratio':>9s}  {'Status':>8s}"
            )
            print(f"  {'-' * 78}")

        rng = np.random.default_rng(42)

        for bid in sorted(reader2._bodies.keys()):
            if not reader1.has_body(bid):
                if verbose:
                    name = BODY_NAMES.get(bid, f"Body {bid}")
                    print(f"  {name:<16s}  {'missing from LEB1':>36s}  {'FAIL':>8s}")
                all_pass = False
                continue

            body = reader2._bodies[bid]
            reference_body = reader1._bodies[bid]
            metadata_errors = []
            if body.coord_type != reference_body.coord_type:
                metadata_errors.append(
                    f"coord_type {body.coord_type} != {reference_body.coord_type}"
                )
            if body.degree != reference_body.degree:
                metadata_errors.append(
                    f"degree {body.degree} != {reference_body.degree}"
                )
            if body.segment_count != reference_body.segment_count:
                metadata_errors.append(
                    "segment_count "
                    f"{body.segment_count} != {reference_body.segment_count}"
                )
            if body.interval_days != reference_body.interval_days:
                metadata_errors.append(
                    "interval_days "
                    f"{body.interval_days} != {reference_body.interval_days}"
                )
            if body.components != reference_body.components:
                metadata_errors.append(
                    f"components {body.components} != {reference_body.components}"
                )
            coverage_values = (
                body.jd_start,
                body.jd_end,
                reference_body.jd_start,
                reference_body.jd_end,
            )
            if not np.all(np.isfinite(coverage_values)):
                metadata_errors.append("non-finite coverage metadata")
            elif expected_group and (
                body.jd_start != reference_body.jd_start
                or body.jd_end != reference_body.jd_end
            ):
                metadata_errors.append(
                    "declared-group coverage "
                    f"[{body.jd_start}, {body.jd_end}] != LEB1 "
                    f"[{reference_body.jd_start}, {reference_body.jd_end}]"
                )
            elif (
                body.jd_start < reference_body.jd_start
                or body.jd_end > reference_body.jd_end
                or body.jd_start >= body.jd_end
            ):
                metadata_errors.append(
                    "coverage "
                    f"[{body.jd_start}, {body.jd_end}] is not contained in LEB1 "
                    f"[{reference_body.jd_start}, {reference_body.jd_end}]"
                )
            if metadata_errors:
                if verbose:
                    name = BODY_NAMES.get(bid, f"Body {bid}")
                    detail = "; ".join(metadata_errors)
                    print(f"  {name:<16s}  metadata mismatch: {detail}  {'FAIL':>8s}")
                all_pass = False
                continue

            jds = rng.uniform(body.jd_start + 0.01, body.jd_end - 0.01, n_samples)

            max_errors = [0.0, 0.0, 0.0]
            component_units = ("?", "?", "?")
            for jd in jds:
                try:
                    pos1, _ = reader1.eval_body(bid, float(jd))
                    pos2, _ = reader2.eval_body(bid, float(jd))
                except Exception:
                    max_errors = [float("inf")] * 3
                    break
                if len(pos1) != 3 or len(pos2) != 3:
                    max_errors = [float("inf")] * 3
                    break
                if not np.all(np.isfinite((*pos1, *pos2))):
                    max_errors = [float("inf")] * 3
                    break
                try:
                    errors, component_units = _native_component_errors(
                        body.coord_type, pos1, pos2
                    )
                except ValueError:
                    max_errors = [float("inf")] * 3
                    break
                if not np.all(np.isfinite(errors)):
                    max_errors = [float("inf")] * 3
                    break
                max_errors = [
                    max(previous, current)
                    for previous, current in zip(max_errors, errors)
                ]

            tolerance = _verification_tolerance(bid, body.degree)
            ratios = [error / tolerance for error in max_errors]
            worst_component = max(range(3), key=ratios.__getitem__)
            max_err = max_errors[worst_component]
            ratio = ratios[worst_component]
            unit = component_units[worst_component]
            passed = all(error <= tolerance for error in max_errors)
            status = "PASS" if passed else "FAIL"
            if not passed:
                all_pass = False

            if verbose:
                name = BODY_NAMES.get(bid, f"Body {bid}")
                print(
                    f"  {name:<16s}  {max_err:>18.2e}  "
                    f"{unit:>5s}  {tolerance:>10.2e}  {ratio:>9.2f}  {status:>8s}"
                )

        reader1.close()
    else:
        if verbose:
            print(f"Verifying {leb2_path} (reader smoke test)")

        jd_mid = (reader2.jd_range[0] + reader2.jd_range[1]) / 2
        for bid in sorted(reader2._bodies.keys()):
            try:
                pos, vel = reader2.eval_body(bid, jd_mid)
                if len(pos) != 3 or len(vel) != 3:
                    raise ValueError("expected three position and velocity components")
                if not np.all(np.isfinite((*pos, *vel))):
                    raise ValueError("non-finite position or velocity")
            except Exception as e:
                if verbose:
                    name = BODY_NAMES.get(bid, f"Body {bid}")
                    print(f"  {name}: FAIL - {e}")
                all_pass = False

        if verbose and all_pass:
            print(f"  All {len(reader2._bodies)} bodies readable: PASS")

    reader2.close()

    if verbose:
        print(f"\n  Result: {'ALL PASS' if all_pass else 'FAILURES DETECTED'}")

    return all_pass


# =============================================================================
# BATCH CONVERSION (all groups for a tier)
# =============================================================================


def convert_all_groups(
    input_path: str,
    output_dir: str,
    tier_name: str,
    target_precision: float = DEFAULT_TARGET_AU,
    verbose: bool = True,
) -> None:
    """Convert a LEB1 file into separate LEB2 files for each group.

    Produces one ``{output_dir}/{tier_name}_{group}.leb2`` per registered
    group, all converted from the merged main ``input_path``.

    Args:
        input_path: Source LEB1 file.
        output_dir: Output directory.
        tier_name: Tier name for file naming (e.g. "base", "medium").
        target_precision: Numeric coefficient target in native stored units.
        verbose: Print progress.
    """
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print(f"=== Batch conversion: {input_path} -> {output_dir}/ ===\n")

    total_size = 0
    converted = 0

    for group_name in LEB2_GROUPS:
        output_path = os.path.join(output_dir, f"{tier_name}_{group_name}.leb2")

        if verbose:
            print(f"\n--- Group: {group_name} ---")

        convert_leb1_to_leb2(
            input_path=input_path,
            output_path=output_path,
            group=group_name,
            expected_tier=tier_name,
            target_precision=target_precision,
            verbose=verbose,
        )

        total_size += os.path.getsize(output_path)
        converted += 1

    input_size = os.path.getsize(input_path)
    if verbose:
        print("\n=== Summary ===")
        print(f"  Input:  {input_size / 1e6:.1f} MB ({os.path.basename(input_path)})")
        print(f"  Output: {total_size / 1e6:.1f} MB (total, {converted} files)")
        if total_size:
            print(f"  Ratio:  {input_size / total_size:.1f}x")


# =============================================================================
# CLI
# =============================================================================


def main():
    """Dispatch LEB2 generation, conversion, verification, and inspection.

    Chunking, compression, body-group boundaries, and CLI defaults are
    LibEphemeris format choices.  Astronomical samples retain the public source
    channels registered by the generator, so conversion does not introduce a
    second scientific authority.
    """
    parser = argparse.ArgumentParser(
        description="Generate or convert LEB2 compressed ephemeris files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Body groups:\n"
            "  core       : Sun-Pluto, Earth, Mean/True Node, Mean Apogee (14 bodies)\n"
            "  asteroids  : Chiron, Ceres, Pallas, Juno, Vesta (5 bodies)\n"
            "  exotics    : Centaurs, TNOs, and NEAs from the canonical registry\n"
            "  apogee     : Oscu Apogee, Interp Apogee/Perigee (3 bodies)\n"
            "\nExamples:\n"
            "  # Convert base tier, core group only:\n"
            "  python generate_leb2.py convert data/leb/ephemeris_base.leb -o core_base.leb --group core --tier base\n"
            "\n"
            "  # Convert all groups from base tier:\n"
            "  python generate_leb2.py convert-all data/leb/ephemeris_base.leb -o data/leb2/ --tier-name base\n"
            "\n"
            "  # Generate from scratch:\n"
            "  python generate_leb2.py generate --tier base --group core -o core_base.leb\n"
            "\n"
            "  # Verify against LEB1 reference:\n"
            "  python generate_leb2.py verify core_base.leb --reference data/leb/ephemeris_base.leb --group core --tier base\n"
        ),
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # --- convert ---
    p_conv = subparsers.add_parser("convert", help="Convert LEB1 -> LEB2")
    p_conv.add_argument("input", help="Input LEB1 file")
    p_conv.add_argument("-o", "--output", required=True, help="Output LEB2 file")
    p_conv.add_argument(
        "--group",
        choices=list(LEB2_GROUPS.keys()),
        help="Only include bodies from this group",
    )
    p_conv.add_argument(
        "--tier",
        choices=["base", "medium", "extended"],
        help="Authenticate tier-specific body inventories",
    )
    p_conv.add_argument(
        "--precision",
        type=float,
        default=DEFAULT_TARGET_AU,
        help=(
            "Numeric coefficient target in native stored units "
            f"(default: {DEFAULT_TARGET_AU})"
        ),
    )
    p_conv.add_argument("-q", "--quiet", action="store_true")

    # --- convert-all ---
    p_all = subparsers.add_parser(
        "convert-all", help="Convert LEB1 -> LEB2 for all groups"
    )
    p_all.add_argument("input", help="Input LEB1 file")
    p_all.add_argument("-o", "--output-dir", required=True, help="Output directory")
    p_all.add_argument(
        "--tier-name",
        choices=["base", "medium", "extended"],
        required=True,
        help="Tier name for file naming (base/medium/extended)",
    )
    p_all.add_argument(
        "--precision",
        type=float,
        default=DEFAULT_TARGET_AU,
        help=(
            "Numeric coefficient target in native stored units "
            f"(default: {DEFAULT_TARGET_AU})"
        ),
    )
    p_all.add_argument("-q", "--quiet", action="store_true")

    # --- generate ---
    p_gen = subparsers.add_parser(
        "generate", help="Generate LEB2 from scratch via Skyfield"
    )
    p_gen.add_argument("-o", "--output", required=True, help="Output LEB2 file")
    p_gen.add_argument(
        "--tier",
        "-t",
        choices=["base", "medium", "extended"],
        required=True,
        help="Precision tier",
    )
    p_gen.add_argument(
        "--group",
        choices=list(LEB2_GROUPS.keys()),
        help="Only include bodies from this group",
    )
    p_gen.add_argument("--start", type=int, help="Override start year")
    p_gen.add_argument("--end", type=int, help="Override end year")
    p_gen.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p_gen.add_argument(
        "--precision",
        type=float,
        default=DEFAULT_TARGET_AU,
        help=(
            "Numeric coefficient target in native stored units "
            f"(default: {DEFAULT_TARGET_AU})"
        ),
    )
    p_gen.add_argument("-q", "--quiet", action="store_true")

    # --- verify ---
    p_ver = subparsers.add_parser("verify", help="Verify a LEB2 file")
    p_ver.add_argument("input", help="LEB2 file to verify")
    p_ver.add_argument(
        "--reference",
        help="LEB1 reference file for comparison",
    )
    p_ver.add_argument(
        "--group",
        choices=list(LEB2_GROUPS),
        help="Require the exact canonical inventory for this named group",
    )
    p_ver.add_argument(
        "--tier",
        choices=["base", "medium", "extended"],
        help="Authenticate tier-specific body inventories",
    )
    p_ver.add_argument("--samples", type=int, default=200, help="Samples per body")
    p_ver.add_argument("-q", "--quiet", action="store_true")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == "convert":
        convert_leb1_to_leb2(
            input_path=args.input,
            output_path=args.output,
            group=args.group,
            expected_tier=args.tier,
            target_precision=args.precision,
            verbose=not args.quiet,
        )

    elif args.command == "convert-all":
        convert_all_groups(
            input_path=args.input,
            output_dir=args.output_dir,
            tier_name=args.tier_name,
            target_precision=args.precision,
            verbose=not args.quiet,
        )

    elif args.command == "generate":
        generate_and_compress(
            output_path=args.output,
            tier=args.tier,
            group=args.group,
            start_year=args.start,
            end_year=args.end,
            workers=args.workers,
            target_precision=args.precision,
            verbose=not args.quiet,
        )

    elif args.command == "verify":
        ok = verify_leb2(
            leb2_path=args.input,
            reference_leb1=args.reference,
            n_samples=args.samples,
            verbose=not args.quiet,
            expected_group=args.group,
            expected_tier=args.tier,
        )
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
