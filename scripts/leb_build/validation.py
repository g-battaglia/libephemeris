"""Structural acceptance of build artifacts, separate from scientific verification.

Provenance:
    Project container validation. Position probes use the local readers; accuracy
    comparison remains exclusively in the existing registered verifiers.
"""

from __future__ import annotations

from contextlib import ExitStack
import math
from pathlib import Path
from typing import Any

from .plan import Job, tier_range
from .storage import sha256_file


def _validate_chunk_layout(reader: Any) -> None:
    """Check every chunk index before any decompression allocation."""
    from libephemeris.leb_format import (
        CHUNK_ENTRY_SIZE,
        CHUNK_INDEX_HEADER_SIZE,
        SECTION_COMPRESSED_CHEBYSHEV,
        segment_byte_size,
    )

    section = reader._sections[SECTION_COMPRESSED_CHEBYSHEV]
    for body_id, body in reader._bodies.items():
        chunks = reader._chunk_index[body_id]
        end = body.data_offset + body.compressed_size
        if (
            not chunks
            or body.data_offset < section.offset
            or end > section.offset + section.size
        ):
            raise ValueError(f"Invalid compressed body bounds: {body_id}")
        cursor = (
            body.data_offset + CHUNK_INDEX_HEADER_SIZE + len(chunks) * CHUNK_ENTRY_SIZE
        )
        segments = 0
        segment_bytes = segment_byte_size(body.degree, body.components)
        for chunk in chunks:
            expected_start = body.jd_start + segments * body.interval_days
            expected_end = expected_start + chunk.segment_count * body.interval_days
            if (
                chunk.segment_start != segments
                or chunk.segment_count <= 0
                or chunk.blob_offset != cursor
                or chunk.compressed_size <= 0
                or chunk.uncompressed_size != chunk.segment_count * segment_bytes
                or not math.isclose(
                    chunk.jd_start, expected_start, rel_tol=0.0, abs_tol=1e-7
                )
                or not math.isclose(
                    chunk.jd_end, expected_end, rel_tol=0.0, abs_tol=1e-7
                )
            ):
                raise ValueError(f"Invalid chunk layout: {body_id}")
            segments += chunk.segment_count
            cursor += chunk.compressed_size
        if (
            segments != body.segment_count
            or cursor != end
            or (body.uncompressed_size != segments * segment_bytes)
        ):
            raise ValueError(f"Invalid chunk totals: {body_id}")


def inspect_artifact(path: Path, job: Job) -> dict[str, Any]:
    """Validate inventory, ranges and storage, then return a content attestation.

    This does not grant scientific PASS. The runner records successful verifier
    jobs separately, and completes the build only when all phases have passed.
    """
    from libephemeris.leb_format import (
        HEADER_SIZE,
        LEB2_VERSION,
        SECTION_DIR_SIZE,
        SECTION_NUTATION,
        SECTION_DELTA_T,
        SECTION_STARS,
    )
    from libephemeris.leb_reader import LEBReader
    from libephemeris.leb2_reader import LEB2Reader
    from scripts.generate_leb import _MergeSource

    is_leb2 = job.kind in ("convert", "verify2")
    with ExitStack() as stack:
        if not is_leb2:
            _MergeSource(str(path), stack)
        reader: LEBReader | LEB2Reader
        if is_leb2:
            reader = stack.enter_context(LEB2Reader(str(path)))
        else:
            reader = stack.enter_context(LEBReader(str(path)))
        header = reader._header
        if is_leb2 and header.version != LEB2_VERSION:
            raise ValueError("Expected current chunked LEB2 version")
        if tuple(sorted(reader._bodies)) != tuple(sorted(job.bodies)) or (
            header.body_count != len(job.bodies)
        ):
            raise ValueError(f"Unexpected body inventory for {job.id}")
        if reader.jd_range != tier_range(job.tier):
            raise ValueError(f"Unexpected global coverage for {job.id}")
        directory_end = HEADER_SIZE + header.section_count * SECTION_DIR_SIZE
        size = path.stat().st_size
        spans = []
        for section in reader._sections.values():
            if section.offset < directory_end or section.offset + section.size > size:
                raise ValueError(f"Invalid section bounds for {job.id}")
            if section.size:
                spans.append((section.offset, section.offset + section.size))
        spans.sort()
        if any(b > c for (_, b), (c, _) in zip(spans, spans[1:])):
            raise ValueError(f"Overlapping sections for {job.id}")
        if is_leb2:
            _validate_chunk_layout(reader)
        auxiliary_ids = {SECTION_NUTATION, SECTION_DELTA_T, SECTION_STARS}
        needs_aux = (job.kind == "generate" and job.bodies == (0,)) or (
            job.kind == "merge" or (is_leb2 and job.group == "core")
        )
        if needs_aux and (
            not auxiliary_ids <= reader._sections.keys()
            or not reader.has_nutation()
            or not reader._delta_t_jds
            or not reader._stars
        ):
            raise ValueError(f"Missing auxiliary data for {job.id}")
        if is_leb2 and job.group != "core" and auxiliary_ids & reader._sections.keys():
            raise ValueError(f"Non-core LEB2 contains shared auxiliary data: {job.id}")
        coverage = {}
        for body_id, body in reader._bodies.items():
            values = (body.jd_start, body.jd_end, body.interval_days)
            if (
                not all(math.isfinite(value) for value in values)
                or body.jd_start >= body.jd_end
                or body.interval_days <= 0
                or body.segment_count <= 0
                or not 1 <= body.components <= 6
                or body.degree > 64
                or body.jd_start < header.jd_start
                or body.jd_end > header.jd_end
            ):
                raise ValueError(f"Invalid body coverage for {job.id}")
            for jd in (body.jd_start, (body.jd_start + body.jd_end) / 2, body.jd_end):
                position, velocity = reader.eval_body(body_id, jd)
                if not all(math.isfinite(value) for value in (*position, *velocity)):
                    raise ValueError(f"Non-finite stored state for {job.id}")
            coverage[str(body_id)] = [float(body.jd_start), float(body.jd_end)]
        return {
            "sha256": sha256_file(path),
            "size": size,
            "bodies": list(sorted(reader._bodies)),
            "jd_range": list(reader.jd_range),
            "body_ranges": coverage,
            "format": "LEB2 v2" if is_leb2 else "LEB1",
        }
