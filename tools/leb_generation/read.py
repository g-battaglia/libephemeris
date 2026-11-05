# SPDX-License-Identifier: AGPL-3.0
"""Checked metadata and bounded read-only sources for project LEB containers."""

from pathlib import Path

from . import format as wire
from .writer import BodyData, SectionData


def _directory(stream, size):
    if size < wire.HEADER_SIZE:
        raise ValueError("truncated LEB header")
    stream.seek(0)
    header = wire.read_header(stream.read(wire.HEADER_SIZE))
    valid = ((header.magic == wire.MAGIC and header.version == wire.VERSION) or
             (header.magic == wire.LEB2_MAGIC and header.version in (wire.LEB2_VERSION_V1, wire.LEB2_VERSION)))
    if not valid:
        raise ValueError("unsupported LEB magic/version")
    wire.validate_entry_count(header.section_count, wire.SECTION_DIR_SIZE,
                              size - wire.HEADER_SIZE, "section directory")
    sections = {}
    for ordinal in range(header.section_count):
        stream.seek(wire.HEADER_SIZE + ordinal * wire.SECTION_DIR_SIZE)
        section = wire.read_section_dir(stream.read(wire.SECTION_DIR_SIZE), 0)
        if section.section_id in sections or section.offset > size or section.size > size - section.offset:
            raise ValueError("duplicate or out-of-bounds source section")
        sections[section.section_id] = section
    coefficient_id = wire.SECTION_CHEBYSHEV if header.magic == wire.MAGIC else wire.SECTION_COMPRESSED_CHEBYSHEV
    if wire.SECTION_BODY_INDEX not in sections or coefficient_id not in sections:
        raise ValueError("source is missing its body index or coefficient section")
    return header, sections, coefficient_id


def _entries(stream, header, sections, coefficient_id):
    compressed = header.magic == wire.LEB2_MAGIC
    width = wire.COMPRESSED_BODY_ENTRY_SIZE if compressed else wire.BODY_ENTRY_SIZE
    decode = wire.read_compressed_body_entry if compressed else wire.read_body_entry
    index, coefficients = sections[wire.SECTION_BODY_INDEX], sections[coefficient_id]
    wire.validate_entry_count(header.body_count, width, index.size, "body index")
    entries = {}
    for ordinal in range(header.body_count):
        stream.seek(index.offset + ordinal * width)
        entry = decode(stream.read(width), 0)
        length = entry.compressed_size if compressed else entry.segment_count * wire.segment_byte_size(entry.degree, entry.components)
        if entry.data_offset < coefficients.offset or entry.data_offset + length > coefficients.offset + coefficients.size:
            raise ValueError(f"body {entry.body_id} coefficients exceed their source section")
        entries[entry.body_id] = entry
    return list(entries.values())


def read_body_index(source: Path):
    """Read only the selected file's LEB1/LEB2 body metadata, never companions.

    Directory/index/payload bounds are checked; actual coefficient and compressed
    chunk validity remain the native reader's responsibility. Repeated IDs keep
    their final physical entry, as in the native index.
    """
    with Path(source).open("rb") as stream:
        stream.seek(0, 2)
        header, sections, coefficient_id = _directory(stream, stream.tell())
        return header, _entries(stream, header, sections, coefficient_id)


def payload_probes(source: Path, header, entries):
    """Select an evaluation inside every stored compressed chunk.

    Warming file pages does not decompress data. Native evaluation at one
    segment midpoint per chunk exercises its real decoder. Chunk selection
    uses segment indices, never the advisory dates in the chunk directory.
    Raw LEB1 and monolithic LEB2 need one probe per body; their payload bounds
    have already been checked. Native admission validates the chunk metadata.
    """
    probes = []
    with Path(source).open("rb") as stream:
        for entry in entries:
            starts = [0]
            if header.magic == wire.LEB2_MAGIC and header.version == wire.LEB2_VERSION:
                stream.seek(entry.data_offset)
                count, _ = wire.read_chunk_index_header(stream.read(wire.CHUNK_INDEX_HEADER_SIZE), 0)
                wire.validate_entry_count(count, wire.CHUNK_ENTRY_SIZE,
                                          entry.compressed_size - wire.CHUNK_INDEX_HEADER_SIZE, "chunk index")
                starts = [wire.read_chunk_entry(stream.read(wire.CHUNK_ENTRY_SIZE), 0).segment_start
                          for _ in range(count)]
            probes.append([entry.body_id, [min(entry.jd_end, entry.jd_start + (start + 0.5) * entry.interval_days)
                                           for start in starts]])
    return probes


def read_leb1(source: Path) -> tuple[wire.FileHeader, list[BodyData], dict[int, SectionData]]:
    """Read indices only; coefficient and auxiliary data remain file-backed.

    Repeated body IDs resolve to their final physical row, matching Rust.
    Section and payload bounds are checked before a caller reads any values.
    No source file is modified, mapped writable or copied into memory.
    """
    source = Path(source)
    with source.open("rb") as stream:
        stream.seek(0, 2)
        size = stream.tell()
        header, sections, coefficient_id = _directory(stream, size)
        if header.magic != wire.MAGIC or header.version != wire.VERSION:
            raise ValueError("operation requires a LEB1 version-1 source")
        entries = _entries(stream, header, sections, coefficient_id)
    auxiliary = {key: SectionData(source, value.offset, value.size) for key, value in sections.items()
                 if key not in (wire.SECTION_BODY_INDEX, wire.SECTION_CHEBYSHEV)}
    return header, [BodyData(entry, path=source) for entry in entries], auxiliary
