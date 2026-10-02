# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Small synthetic native artifacts shared by unit and PostgreSQL tests."""

from __future__ import annotations

import struct
from pathlib import Path

from libephemeris import leb_format as layout
from libephemeris.leb_compression import compress_body
from scripts.generate_leb2 import write_leb2, write_leb2_chunked


def make_artifacts(directory: Path) -> list[Path]:
    """Create LEB1, LEB2 v1 and LEB2 v2 with identical synthetic coefficients.

    Args:
        directory: Disposable directory owned by the test.

    Returns:
        Three artifact paths containing body, nutation, Delta-T and star data.
    """
    body = layout.BodyEntry(0, 0, 2, 10.0, 14.0, 2.0, 2, 3, 0)
    coefficients = struct.pack(
        "<18d", *([1.0, 2.0, 3.0, -4.0, 0.5, 0.0, 7.0, 0.0, 0.0] * 2)
    )
    nutation = struct.pack(layout.NUTATION_HEADER_FMT, 10.0, 14.0, 4.0, 2, 2, 1, 0)
    nutation += struct.pack("<6d", 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    delta_t = struct.pack("<II4d", 2, 0, 10.0, 0.1, 14.0, 0.2)
    stars = struct.pack(
        layout.STAR_ENTRY_FMT, 42, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, b"\x00" * 4
    )

    section_ids = [
        layout.SECTION_BODY_INDEX,
        layout.SECTION_CHEBYSHEV,
        layout.SECTION_NUTATION,
        layout.SECTION_DELTA_T,
        layout.SECTION_STARS,
    ]
    directory_end = layout.HEADER_SIZE + len(section_ids) * layout.SECTION_DIR_SIZE
    body.data_offset = directory_end + layout.BODY_ENTRY_SIZE
    index = bytearray(layout.BODY_ENTRY_SIZE)
    layout.write_body_entry(index, 0, body)
    payloads = [bytes(index), coefficients, nutation, delta_t, stars]
    buffer = bytearray(directory_end + sum(len(payload) for payload in payloads))
    layout.write_header(
        buffer, layout.FileHeader(layout.MAGIC, 1, 5, 1, 10.0, 14.0, 0.0, 0)
    )
    offset = directory_end
    for index, (section_id, payload) in enumerate(zip(section_ids, payloads)):
        section = layout.SectionEntry(section_id, offset, len(payload))
        layout.write_section_dir(
            buffer, layout.HEADER_SIZE + index * layout.SECTION_DIR_SIZE, section
        )
        buffer[offset : offset + len(payload)] = payload
        offset += len(payload)
    leb1 = directory / "base_synthetic.leb"
    leb1.write_bytes(buffer)

    # Keeping all 52 mantissa bits ensures these are format tests, not tests
    # of compression approximation. Both codecs reuse the project writer.
    compressed = compress_body(coefficients, 2, 2, 3, [52, 52, 52])
    leb2_v1 = directory / "base_synthetic_v1.leb2"
    write_leb2(
        str(leb2_v1),
        [(body, compressed, len(coefficients))],
        nutation,
        delta_t,
        stars,
        10.0,
        14.0,
        0.0,
        verbose=False,
    )
    chunks = []
    for index in range(2):
        payload = coefficients[index * 72 : (index + 1) * 72]
        chunks.append((compress_body(payload, 1, 2, 3, [52, 52, 52]), len(payload)))
    leb2_v2 = directory / "base_synthetic_v2.leb2"
    write_leb2_chunked(
        str(leb2_v2),
        [(body, chunks, len(coefficients))],
        nutation,
        delta_t,
        stars,
        10.0,
        14.0,
        0.0,
        chunk_interval_days=2.0,
        verbose=False,
    )
    return [leb1, leb2_v1, leb2_v2]
