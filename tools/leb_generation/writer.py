# SPDX-License-Identifier: AGPL-3.0
"""Bounded-memory, atomic LEB2 production for offline generation/conversion."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import compression
from . import format as wire


@dataclass(frozen=True)
class BodyData:
    """One body's metadata and explicit binary64 coefficient source.

    File coefficients begin at entry.data_offset and are little-endian,
    segment/component/order major. An in-memory array has the same axes.
    Exactly one source is required; neither source invokes an ephemeris engine.
    """

    entry: wire.BodyEntry
    path: Path | None = None
    coefficients: np.ndarray | None = None

    def chunks(self, count: int) -> Iterator[tuple[int, np.ndarray]]:
        entry = self.entry
        shape = (entry.segment_count, entry.components, entry.degree + 1)
        if count < 1 or min(shape) < 1 or not math.isfinite(entry.interval_days) or entry.interval_days <= 0:
            raise ValueError("invalid coefficient shape or temporal chunk interval")
        if (self.path is None) == (self.coefficients is None):
            raise ValueError("exactly one coefficient source is required")
        if self.coefficients is not None:
            values = self.coefficients
            if values.shape != shape or values.dtype.kind != "f" or values.dtype.itemsize != 8:
                raise ValueError("coefficients must have the declared binary64 shape")
            for start in range(0, shape[0], count):
                yield start, np.asarray(values[start:start + count], dtype="<f8")
            return
        segment_bytes = wire.segment_byte_size(entry.degree, entry.components)
        with self.path.open("rb") as source:
            source.seek(entry.data_offset)
            for start in range(0, shape[0], count):
                length = min(count, shape[0] - start)
                raw = source.read(length * segment_bytes)
                if len(raw) != length * segment_bytes:
                    raise ValueError(f"truncated coefficient source for body {entry.body_id}")
                yield start, np.frombuffer(raw, dtype="<f8").reshape(length, *shape[1:])


def compressed_chunk(coefficients: np.ndarray, body_id: int, *, lossless: bool, target: float | None) -> bytes:
    """Compress one chunk and verify its native-component polynomial bound.

    The sum bound follows |T_k(t)| <= 1 on [-1,1]. It is not a guarantee on
    transformed apparent coordinates; complete generation must additionally
    verify those coordinates against the actual native source provider.
    """
    if not np.isfinite(coefficients).all():
        raise ValueError(f"nonfinite coefficients for body {body_id}")
    degree = coefficients.shape[2] - 1
    budget = compression.BODY_TARGET_AU.get(body_id, compression.DEFAULT_TARGET_AU)
    if target is not None:
        if not math.isfinite(target) or target <= 0:
            raise ValueError("coefficient precision must be finite and positive")
        budget = min(budget, target)
    if not math.isfinite(budget) or budget <= 0:
        raise ValueError("coefficient precision must be finite and positive")
    bits = [52] * (degree + 1) if lossless else compression.compute_mantissa_bits(coefficients, budget / (degree + 1))
    raw = coefficients.tobytes()
    blob = compression.compress_body(raw, coefficients.shape[0], degree, coefficients.shape[1], bits)
    decoded = compression.decompress_body(blob, len(raw), coefficients.shape[0], degree, coefficients.shape[1])
    if lossless:
        if decoded != raw:
            raise ValueError("lossless coefficient verification failed")
    else:
        restored = np.frombuffer(decoded, dtype="<f8").reshape(coefficients.shape)
        bound = np.sum(np.abs(restored - coefficients), axis=2)
        if np.max(bound) > budget:
            raise ValueError(f"coefficient sum error exceeds the fixed budget for body {body_id}")
    return blob


def verify_native(path: Path, bodies: list[BodyData]) -> None:
    """Admit the complete container and force every native compressed chunk read."""
    verify_native_index(path, [body.entry for body in bodies])


def verify_native_index(path: Path, entries) -> None:
    """Verify the selected container's actual body index and every payload chunk."""
    from .read import payload_probes, read_body_index

    header, actual_entries = read_body_index(path)
    # A fresh process selects the indexed Rust reader before package import.
    # This preserves the caller's runtime and bounds native validation memory,
    # including for multi-gigabyte output and all compressed temporal chunks.
    program = """
import json, sys
import libephemeris as e
import libephemeris._core
e.set_leb_file(sys.argv[1])
reader = e.get_leb_reader()
expected, probes = json.load(sys.stdin)
assert reader is not None
assert set(reader.body_ids()) == set(expected)
for body, epochs in probes:
    for epoch in epochs:
        reader.eval_body(body, epoch)
        reader.cool()
reader.close()
"""
    with tempfile.TemporaryDirectory(prefix="libephemeris-leb-check-") as directory:
        root = Path(directory)
        environment = {key: value for key, value in os.environ.items() if not key.startswith("LIBEPHEMERIS_")}
        environment.update(
            LIBEPHEMERIS_MODE="leb", LIBEPHEMERIS_LEB_READER="indexed",
            LIBEPHEMERIS_DATA_DIR=str(root), LIBEPHEMERIS_NETWORK_POLICY="sealed",
            LIBEPHEMERIS_CONFIG=str(root / "absent.toml"), LIBEPHEMERIS_ENV_FILE=str(root / "absent.env"),
            LIBEPHEMERIS_IERS_AUTO_DOWNLOAD="0", LIBEPHEMERIS_AUTO_SPK="0",
        )
        result = subprocess.run(
            [sys.executable, "-I", "-B", "-c", program, str(path.resolve())],
            input=json.dumps([[entry.body_id for entry in entries], payload_probes(path, header, actual_entries)]),
            text=True, capture_output=True, env=environment, check=False,
        )
        if result.returncode:
            raise ValueError(f"native generated-file verification failed: {result.stderr.strip()}")


@dataclass(frozen=True)
class SectionData:
    """An unchanged auxiliary section copied from a read-only source in blocks."""

    path: Path
    offset: int
    size: int

    def copy_to(self, destination) -> None:
        with self.path.open("rb") as source:
            source.seek(self.offset)
            remaining = self.size
            while remaining:
                block = source.read(min(remaining, 1024 * 1024))
                if not block:
                    raise ValueError("truncated auxiliary source section")
                destination.write(block)
                remaining -= len(block)


def write_leb2(
    output: Path,
    bodies: list[BodyData],
    header: wire.FileHeader,
    auxiliary: dict[int, bytes | SectionData] | None = None,
    *,
    chunk_interval_days: float = wire.CHUNK_INTERVAL_DAYS,
    lossless: bool = True,
    target: float | None = None,
    overwrite: bool = False,
    validator=None,
) -> None:
    """Stream LEB2 v2, verify through Rust, then atomically publish the file.

    Input and output paths must differ. Existing output is preserved unless
    overwrite is selected; any serialization, compression or verification
    failure removes only this operation's temporary file. Chunk memory is
    bounded by chunk_interval_days and the source body's segment interval.
    Auxiliary sections are preserved as opaque project-format bytes.
    """
    output = Path(output)
    if not bodies or len({body.entry.body_id for body in bodies}) != len(bodies):
        raise ValueError("generation requires a nonempty unique body inventory")
    if not math.isfinite(chunk_interval_days) or chunk_interval_days <= 0:
        raise ValueError("chunk interval must be finite and positive")
    if any(body.path is not None and body.path.resolve() == output.resolve() for body in bodies):
        raise ValueError("output must differ from every coefficient source")
    if output.exists() and not overwrite:
        raise FileExistsError(output)
    auxiliary = dict(auxiliary or {})
    if any(isinstance(data, SectionData) and data.path.resolve() == output.resolve() for data in auxiliary.values()):
        raise ValueError("output must differ from every auxiliary source")
    if set(auxiliary) & {wire.SECTION_BODY_INDEX, wire.SECTION_CHEBYSHEV, wire.SECTION_COMPRESSED_CHEBYSHEV}:
        raise ValueError("auxiliary sections cannot replace body or coefficient indices")
    section_count = 2 + len(auxiliary)
    index_at = wire.HEADER_SIZE + section_count * wire.SECTION_DIR_SIZE
    coefficients_at = index_at + len(bodies) * wire.COMPRESSED_BODY_ENTRY_SIZE
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{output.name}-", suffix=".tmp", dir=output.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w+b") as destination:
            destination.write(bytes(coefficients_at))
            entries = []
            for body in bodies:
                entry = body.entry
                if not math.isfinite(entry.interval_days) or entry.interval_days <= 0:
                    raise ValueError("body interval must be finite and positive")
                count = max(1, math.floor(chunk_interval_days / entry.interval_days))
                chunks = (entry.segment_count + count - 1) // count
                body_at = destination.tell()
                directory = bytearray(wire.CHUNK_INDEX_HEADER_SIZE + chunks * wire.CHUNK_ENTRY_SIZE)
                wire.write_chunk_index_header(directory, 0, chunks, chunk_interval_days)
                destination.write(directory)
                uncompressed = 0
                for ordinal, (start, values) in enumerate(body.chunks(count)):
                    blob = compressed_chunk(values, entry.body_id, lossless=lossless, target=target)
                    blob_at = destination.tell()
                    destination.write(blob)
                    uncompressed += values.nbytes
                    chunk = wire.ChunkEntry(
                        entry.jd_start + start * entry.interval_days,
                        entry.jd_start + (start + len(values)) * entry.interval_days,
                        blob_at, len(blob), values.nbytes, start, len(values),
                    )
                    wire.write_chunk_entry(directory, wire.CHUNK_INDEX_HEADER_SIZE + ordinal * wire.CHUNK_ENTRY_SIZE, chunk)
                end = destination.tell()
                destination.seek(body_at)
                destination.write(directory)
                destination.seek(end)
                entries.append(wire.CompressedBodyEntry(
                    entry.body_id, entry.coord_type, entry.segment_count,
                    entry.jd_start, entry.jd_end, entry.interval_days, entry.degree,
                    entry.components, body_at, end - body_at, uncompressed,
                ))
            sections = [
                wire.SectionEntry(wire.SECTION_BODY_INDEX, index_at, len(entries) * wire.COMPRESSED_BODY_ENTRY_SIZE),
                wire.SectionEntry(wire.SECTION_COMPRESSED_CHEBYSHEV, coefficients_at, destination.tell() - coefficients_at),
            ]
            for section_id, data in sorted(auxiliary.items()):
                size = data.size if isinstance(data, SectionData) else len(data)
                sections.append(wire.SectionEntry(section_id, destination.tell(), size))
                if isinstance(data, SectionData):
                    data.copy_to(destination)
                else:
                    destination.write(data)
            prefix = bytearray(coefficients_at)
            wire.write_header(prefix, wire.FileHeader(
                wire.LEB2_MAGIC, wire.LEB2_VERSION, section_count, len(bodies),
                header.jd_start, header.jd_end, header.generation_epoch,
                wire.COMPRESSION_ZSTD_TRUNC_SHUFFLE,
            ))
            for ordinal, section in enumerate(sections):
                wire.write_section_dir(prefix, wire.HEADER_SIZE + ordinal * wire.SECTION_DIR_SIZE, section)
            for ordinal, entry in enumerate(entries):
                wire.write_compressed_body_entry(prefix, index_at + ordinal * wire.COMPRESSED_BODY_ENTRY_SIZE, entry)
            destination.seek(0)
            destination.write(prefix)
            destination.flush()
            os.fsync(destination.fileno())
        verify_native(temporary, bodies)
        if validator is not None:
            validator(temporary)
        if overwrite:
            os.replace(temporary, output)
        else:
            # An atomic link also refuses an output created by another writer
            # after the initial existence check. Both files share this directory.
            os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
