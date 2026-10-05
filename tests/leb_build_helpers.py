"""Synthetic project-native LEB fixtures; no external ephemeris samples."""

from __future__ import annotations

from pathlib import Path
import struct

from libephemeris import leb_format as fmt


def write_synthetic_leb(
    path: Path,
    bodies: tuple[int, ...] = (0,),
    *,
    aux: bool = True,
    start: float = 100.0,
    end: float = 110.0,
    segments: int = 2,
) -> None:
    """Write deterministic toy polynomials and optional auxiliary records."""
    body_bytes = segments * fmt.segment_byte_size(1, 3)
    nutation_size = fmt.NUTATION_HEADER_SIZE + (32 if aux else 0)
    sizes = (
        len(bodies) * fmt.BODY_ENTRY_SIZE,
        len(bodies) * body_bytes,
        nutation_size,
        fmt.DELTA_T_HEADER_SIZE + (32 if aux else 0),
        fmt.STAR_ENTRY_SIZE if aux else 0,
    )
    offset = fmt.HEADER_SIZE + len(sizes) * fmt.SECTION_DIR_SIZE
    sections = []
    for sid, size in enumerate(sizes):
        sections.append(fmt.SectionEntry(sid, offset, size))
        offset += size
    data = bytearray(offset)
    fmt.write_header(
        data,
        fmt.FileHeader(
            fmt.MAGIC,
            fmt.VERSION,
            len(sizes),
            len(bodies),
            start,
            end,
            start,
            0,
        ),
    )
    for sid, section in enumerate(sections):
        fmt.write_section_dir(
            data, fmt.HEADER_SIZE + sid * fmt.SECTION_DIR_SIZE, section
        )
    for i, body in enumerate(bodies):
        payload_offset = sections[1].offset + i * body_bytes
        fmt.write_body_entry(
            data,
            sections[0].offset + i * fmt.BODY_ENTRY_SIZE,
            fmt.BodyEntry(
                body,
                fmt.COORD_ICRS_BARY,
                segments,
                start,
                end,
                (end - start) / segments,
                1,
                3,
                payload_offset,
            ),
        )
        for segment in range(segments):
            struct.pack_into(
                "<6d",
                data,
                payload_offset + segment * 48,
                0.1 + body,
                0.001,
                0.2,
                0.002,
                0.3,
                0.003,
            )
    fmt.write_nutation_header(
        data,
        sections[2].offset,
        fmt.NutationHeader(
            start,
            end,
            end - start,
            1,
            2,
            int(aux),
            0,
        ),
    )
    if aux:
        struct.pack_into(
            "<4d",
            data,
            sections[2].offset + fmt.NUTATION_HEADER_SIZE,
            0.01,
            0.001,
            0.02,
            0.002,
        )
    struct.pack_into("<II", data, sections[3].offset, 2 if aux else 0, 0)
    if aux:
        struct.pack_into(
            "<4d",
            data,
            sections[3].offset + fmt.DELTA_T_HEADER_SIZE,
            start,
            0.001,
            end,
            0.002,
        )
        fmt.write_star_entry(
            data,
            sections[4].offset,
            fmt.StarEntry(1, 10.0, 20.0, 0.0, 0.0, 0.1, 0.0, 1.0),
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
