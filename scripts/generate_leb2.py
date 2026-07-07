#!/usr/bin/env python3
"""
LEB2 compressed ephemeris file generator.

Two modes of operation:

1. Convert: Read an existing LEB1 file and produce a compressed LEB2 file.
   python scripts/generate_leb2.py convert data/leb/ephemeris_base.leb -o core_base.leb --group core

2. Generate: Generate from scratch via Skyfield + Chebyshev fitting, then compress.
   python scripts/generate_leb2.py generate --tier base --group core -o core_base.leb

The LEB2 format uses error-bounded lossy compression (mantissa truncation +
coefficient-major reorder + byte shuffle + zstd) achieving 5-15x compression
per body while maintaining <0.001" precision.

Body groups:
  core       : Sun, Moon, Mercury-Pluto, Earth, Mean/True Node, Mean Apogee (14 bodies)
  asteroids  : Chiron, Ceres, Pallas, Juno, Vesta (5 bodies)
  apogee     : Osculating Apogee, Interpolated Apogee/Perigee (3 bodies)
  uranians   : Cupido-Transpluto (9 bodies)
"""

from __future__ import annotations

import argparse
import mmap
import os
import sys
import time
from typing import Optional

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
from libephemeris.exotic_bodies import (
    EXOTIC_IDS,
    name_map as _exotic_names,
)

# =============================================================================
# BODY GROUPS
# =============================================================================

LEB2_GROUPS = {
    "core": sorted([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 10, 11, 12]),
    "asteroids": sorted([15, 17, 18, 19, 20]),
    "apogee": sorted([13, 21, 22]),
    "uranians": sorted([40, 41, 42, 43, 44, 45, 46, 47, 48]),
    "exotics": list(EXOTIC_IDS),  # centaurs/TNOs/NEAs — registry source of truth
}

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
        sec_entries.append(SectionEntry(SECTION_NUTATION, nut_offset, len(nutation_data)))
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
            buf[blob_cursor: blob_cursor + len(blob)] = blob
            blob_cursor += len(blob)
            seg_cursor += chunk_seg_count

        data_cursor += body_data_size

    # Write auxiliary sections
    if nutation_data:
        buf[nut_offset: nut_offset + len(nutation_data)] = nutation_data
    if delta_t_data:
        buf[dt_offset: dt_offset + len(delta_t_data)] = delta_t_data
    if star_data:
        buf[star_offset: star_offset + len(star_data)] = star_data

    # Flush
    os.makedirs(os.path.dirname(os.path.abspath(output)) or ".", exist_ok=True)
    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print(f"\n  Output: {output}")
        print(f"  File size: {total_size / 1e6:.2f} MB ({total_size:,} bytes)")
        print(f"  Format: LEB2 v2 (chunked, {chunk_interval_days / 365.25:.0f}-year chunks)")


# =============================================================================
# CONVERT MODE: LEB1 -> LEB2 (v2 chunked)
# =============================================================================


def convert_leb1_to_leb2(
    input_path: str,
    output_path: str,
    group: Optional[str] = None,
    target_precision: float = DEFAULT_TARGET_AU,
    verbose: bool = True,
) -> None:
    """Convert an existing LEB1 file to LEB2 v2 format (chunked).

    Args:
        input_path: Path to source LEB1 file.
        output_path: Path to output LEB2 file.
        group: Optional body group filter (core/asteroids/apogee/uranians).
        target_precision: Precision target in AU for mantissa truncation.
        verbose: Print progress.
    """
    if verbose:
        chunk_years = CHUNK_INTERVAL_DAYS / 365.25
        print(f"Converting {input_path} -> {output_path} (v2, {chunk_years:.0f}-year chunks)")
        if group:
            print(f"  Group: {group} ({len(LEB2_GROUPS[group])} bodies)")

    t0 = time.time()
    src = LEB1Source(input_path)

    # Select bodies
    if group:
        body_ids = [bid for bid in LEB2_GROUPS[group] if bid in src.bodies]
    else:
        body_ids = sorted(src.bodies.keys())

    if verbose:
        print(f"  Bodies: {len(body_ids)}")
        input_size = os.path.getsize(input_path)
        print(f"  Input size: {input_size / 1e6:.1f} MB")

    # Compress each body into chunks
    body_chunks = []
    total_raw = 0
    total_compressed = 0

    if verbose:
        print(f"\n  {'Body':<16s}  {'Raw KB':>8s}  {'Comp KB':>8s}  {'Ratio':>6s}  {'Chunks':>6s}")
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
            raw, entry.segment_count, entry.degree, entry.components, bits, segs_per_chunk
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

    # Get auxiliary data
    nutation_data = src.get_section_data(SECTION_NUTATION)
    delta_t_data = src.get_section_data(SECTION_DELTA_T)
    star_data = src.get_section_data(SECTION_STARS)

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
    "extended": ("de441.bsp", -5000, 5000, "ephemeris_extended"),
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
        target_precision: Precision target in AU.
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


def verify_leb2(
    leb2_path: str,
    reference_leb1: Optional[str] = None,
    n_samples: int = 200,
    verbose: bool = True,
) -> bool:
    """Verify a LEB2 file against a LEB1 reference or Skyfield.

    Args:
        leb2_path: Path to LEB2 file to verify.
        reference_leb1: Optional LEB1 reference file. If provided, compares
            LEB2 eval_body results against LEB1. Otherwise uses Skyfield.
        n_samples: Number of random JDs to sample per body.
        verbose: Print results.

    Returns:
        True if all bodies pass precision checks.
    """
    from libephemeris.leb2_reader import LEB2Reader

    reader2 = LEB2Reader(leb2_path)
    all_pass = True

    if reference_leb1:
        from libephemeris.leb_reader import LEBReader

        reader1 = LEBReader(reference_leb1)

        if verbose:
            print(f"Verifying {leb2_path} against {reference_leb1}")
            print(f"  Samples per body: {n_samples}")
            arcsec_hdr = 'Max err (")'
            print(
                f"\n  {'Body':<16s}  {'Max err (AU)':>14s}  {arcsec_hdr:>12s}  {'Status':>8s}"
            )
            print(f"  {'-' * 56}")

        rng = np.random.default_rng(42)

        for bid in sorted(reader2._bodies.keys()):
            if not reader1.has_body(bid):
                continue

            body = reader2._bodies[bid]
            jds = rng.uniform(body.jd_start + 0.01, body.jd_end - 0.01, n_samples)

            max_err = 0.0
            for jd in jds:
                pos1, _ = reader1.eval_body(bid, float(jd))
                pos2, _ = reader2.eval_body(bid, float(jd))
                err = max(abs(a - b) for a, b in zip(pos1, pos2))
                max_err = max(max_err, err)

            # Convert to arcsec (rough: 1 AU ≈ 206265")
            err_arcsec = max_err * 206265
            passed = err_arcsec < 1.0  # generous threshold
            status = "PASS" if passed else "FAIL"
            if not passed:
                all_pass = False

            if verbose:
                name = BODY_NAMES.get(bid, f"Body {bid}")
                print(
                    f"  {name:<16s}  {max_err:>14.2e}  {err_arcsec:>12.4f}  {status:>8s}"
                )

        reader1.close()
    else:
        if verbose:
            print(f"Verifying {leb2_path} (reader smoke test)")

        jd_mid = (reader2.jd_range[0] + reader2.jd_range[1]) / 2
        for bid in sorted(reader2._bodies.keys()):
            try:
                pos, vel = reader2.eval_body(bid, jd_mid)
                assert len(pos) == 3
                assert len(vel) == 3
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

    Produces:
      {output_dir}/{tier_name}_core.leb2
      {output_dir}/{tier_name}_asteroids.leb2
      {output_dir}/{tier_name}_apogee.leb2
      {output_dir}/{tier_name}_uranians.leb2

    Args:
        input_path: Source LEB1 file.
        output_dir: Output directory.
        tier_name: Tier name for file naming (e.g. "base", "medium").
        target_precision: Precision target in AU.
        verbose: Print progress.
    """
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print(f"=== Batch conversion: {input_path} -> {output_dir}/ ===\n")

    total_size = 0

    for group_name, body_ids in LEB2_GROUPS.items():
        output_path = os.path.join(output_dir, f"{tier_name}_{group_name}.leb2")

        if verbose:
            print(f"\n--- Group: {group_name} ({len(body_ids)} bodies) ---")

        convert_leb1_to_leb2(
            input_path=input_path,
            output_path=output_path,
            group=group_name,
            target_precision=target_precision,
            verbose=verbose,
        )

        total_size += os.path.getsize(output_path)

    input_size = os.path.getsize(input_path)
    if verbose:
        print("\n=== Summary ===")
        print(f"  Input:  {input_size / 1e6:.1f} MB ({os.path.basename(input_path)})")
        print(f"  Output: {total_size / 1e6:.1f} MB (total, {len(LEB2_GROUPS)} files)")
        print(f"  Ratio:  {input_size / total_size:.1f}x")


# =============================================================================
# CLI
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Generate or convert LEB2 compressed ephemeris files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Body groups:\n"
            "  core       : Sun-Pluto, Earth, Mean/True Node, Mean Apogee (14 bodies)\n"
            "  asteroids  : Chiron, Ceres, Pallas, Juno, Vesta (5 bodies)\n"
            "  apogee     : Oscu Apogee, Interp Apogee/Perigee (3 bodies)\n"
            "  uranians   : Cupido-Transpluto (9 bodies)\n"
            "\nExamples:\n"
            "  # Convert base tier, core group only:\n"
            "  python generate_leb2.py convert data/leb/ephemeris_base.leb -o core_base.leb --group core\n"
            "\n"
            "  # Convert all groups from base tier:\n"
            "  python generate_leb2.py convert-all data/leb/ephemeris_base.leb -o data/leb2/ --tier-name base\n"
            "\n"
            "  # Generate from scratch:\n"
            "  python generate_leb2.py generate --tier base --group core -o core_base.leb\n"
            "\n"
            "  # Verify against LEB1 reference:\n"
            "  python generate_leb2.py verify core_base.leb --reference data/leb/ephemeris_base.leb\n"
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
        "--precision",
        type=float,
        default=DEFAULT_TARGET_AU,
        help=f"Precision target in AU (default: {DEFAULT_TARGET_AU})",
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
        required=True,
        help="Tier name for file naming (base/medium/extended)",
    )
    p_all.add_argument(
        "--precision",
        type=float,
        default=DEFAULT_TARGET_AU,
        help=f"Precision target in AU (default: {DEFAULT_TARGET_AU})",
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
        help=f"Precision target in AU (default: {DEFAULT_TARGET_AU})",
    )
    p_gen.add_argument("-q", "--quiet", action="store_true")

    # --- verify ---
    p_ver = subparsers.add_parser("verify", help="Verify a LEB2 file")
    p_ver.add_argument("input", help="LEB2 file to verify")
    p_ver.add_argument(
        "--reference",
        help="LEB1 reference file for comparison",
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
        )
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
