# SPDX-License-Identifier: AGPL-3.0
"""Stream and atomically publish native-readable LEB1 files and merged parts."""

import os
import hashlib
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import format as wire
from .read import read_leb1
from .source import check_input_hashes
from .writer import BodyData, SectionData, verify_native


def write_leb1(output: Path, bodies: list[BodyData], header: wire.FileHeader,
               auxiliary: dict[int, bytes | SectionData] | None = None, *, overwrite=False, validator=None) -> None:
    """Write bounded coefficient blocks; verify before atomic publication.

    Every coefficient must be finite binary64. Existing data and all input
    files survive serialization, validation or publication failures. The same
    project-format auxiliary sections used by LEB2 are accepted unchanged.
    """
    output = Path(output)
    auxiliary = dict(auxiliary or {})
    if not bodies or len({body.entry.body_id for body in bodies}) != len(bodies):
        raise ValueError("generation requires a nonempty unique body inventory")
    if set(auxiliary) & {wire.SECTION_BODY_INDEX, wire.SECTION_CHEBYSHEV, wire.SECTION_COMPRESSED_CHEBYSHEV}:
        raise ValueError("auxiliary sections cannot replace body or coefficient indices")
    sources = [body.path for body in bodies if body.path is not None]
    sources.extend(data.path for data in auxiliary.values() if isinstance(data, SectionData))
    if any(source.resolve() == output.resolve() for source in sources):
        raise ValueError("output must differ from every source")
    if output.exists() and not overwrite:
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    section_count = 2 + len(auxiliary)
    index_at = wire.HEADER_SIZE + section_count * wire.SECTION_DIR_SIZE
    coefficients_at = index_at + len(bodies) * wire.BODY_ENTRY_SIZE
    descriptor, name = tempfile.mkstemp(prefix=f".{output.name}-", suffix=".tmp", dir=output.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w+b") as destination:
            destination.write(bytes(coefficients_at))
            entries = []
            for body in bodies:
                entries.append(replace(body.entry, data_offset=destination.tell()))
                for _, coefficients in body.chunks(256):
                    if not np.isfinite(coefficients).all():
                        raise ValueError(f"nonfinite coefficients for body {body.entry.body_id}")
                    destination.write(coefficients.tobytes())
            sections = [
                wire.SectionEntry(wire.SECTION_BODY_INDEX, index_at, len(entries) * wire.BODY_ENTRY_SIZE),
                wire.SectionEntry(wire.SECTION_CHEBYSHEV, coefficients_at, destination.tell() - coefficients_at),
            ]
            for section_id, data in sorted(auxiliary.items()):
                size = data.size if isinstance(data, SectionData) else len(data)
                sections.append(wire.SectionEntry(section_id, destination.tell(), size))
                if isinstance(data, SectionData):
                    data.copy_to(destination)
                else:
                    destination.write(data)
            prefix = bytearray(coefficients_at)
            wire.write_header(prefix, replace(header, magic=wire.MAGIC, version=wire.VERSION,
                                               section_count=section_count, body_count=len(bodies), flags=0))
            for ordinal, section in enumerate(sections):
                wire.write_section_dir(prefix, wire.HEADER_SIZE + ordinal * wire.SECTION_DIR_SIZE, section)
            for ordinal, entry in enumerate(entries):
                wire.write_body_entry(prefix, index_at + ordinal * wire.BODY_ENTRY_SIZE, entry)
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
            os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def merge_leb1(inputs: list[Path], output: Path, *, overwrite=False) -> None:
    """Merge disjoint body inventories with the same global domain.

    Ambiguous duplicate bodies or differing shared auxiliary bytes are errors;
    there is no implicit winner. Coefficients and sections are streamed from
    the original files without recomputation, truncation or clock changes.
    """
    if not inputs:
        raise ValueError("merge requires at least one input")
    parts = [read_leb1(source) for source in inputs]
    header = parts[0][0]
    bodies, auxiliary = [], {}
    identifiers = set()
    for candidate, rows, sections in parts:
        if (candidate.jd_start, candidate.jd_end) != (header.jd_start, header.jd_end):
            raise ValueError("merged files must have the same global date range")
        for body in rows:
            if body.entry.body_id in identifiers:
                raise ValueError(f"duplicate merge body: {body.entry.body_id}")
            identifiers.add(body.entry.body_id)
            bodies.append(body)
        for section_id, section in sections.items():
            if section_id in auxiliary and _section_digest(section) != _section_digest(auxiliary[section_id]):
                raise ValueError(f"conflicting shared auxiliary section: {section_id}")
            auxiliary[section_id] = section
    write_leb1(output, bodies, header, auxiliary, overwrite=overwrite)


def _check_replacement_paths(main: Path, inputs: list[Path], output: Path, overwrite: bool) -> None:
    """Protect every acquired input name, including symlink and hardlink aliases."""
    def aliases(left, right):
        return (left.resolve() == right.resolve() or
                (left.exists() and right.exists() and left.samefile(right)))

    sources = [main, *inputs]
    for ordinal, source in enumerate(sources):
        if any(aliases(source, previous) for previous in sources[:ordinal]):
            raise ValueError("replacement source files must be distinct")
    if any(aliases(output, source) for source in inputs):
        raise ValueError("replacement output must differ from every partial source")
    if aliases(output, main) and not overwrite:
        raise FileExistsError("replacing the main file requires overwrite=True")
    if os.path.lexists(output) and not overwrite:
        raise FileExistsError(output)


def replace_leb1(main: Path, inputs: list[Path], output: Path, *,
                 expected_bodies: set[int], overwrite=False) -> None:
    """Replace exactly the selected bodies, retaining the main file's auxiliaries.

    All sources must be distinct immutable LEB1 images with identical global
    Julian-date endpoints. Replacement parts must supply the complete nonempty
    expected inventory without duplicates. Their auxiliary sections are unused;
    every main auxiliary section, including opaque payloads, is copied unchanged.
    Unselected body records and coefficient words retain their original values,
    apart from relocated payload offsets; resulting bodies are sorted by ID.

    A sibling staged file passes native validation before atomic publication.
    In-place replacement requires overwrite=True; partial sources are never
    destinations. Source hashes are checked after verification, immediately
    before publication. Concurrent source mutation is unsupported, rather than
    being made into a snapshot by these checks. Failures preserve the destination.
    An existing destination retains its ordinary read/write/execute permissions.
    """
    main, output = Path(main), Path(output)
    inputs = [Path(path) for path in inputs]
    if not inputs or not expected_bodies:
        raise ValueError("replacement requires a nonempty expected body inventory and partials")
    _check_replacement_paths(main, inputs, output, overwrite)
    permissions = output.stat().st_mode & 0o777 if output.exists() else None
    identities = []
    for path in [main, *inputs]:
        with path.open("rb") as source:
            identities.append((path, hashlib.file_digest(source, "sha256").hexdigest()))

    header, original, auxiliary = read_leb1(main)
    if header.body_count != len(original):
        raise ValueError("duplicate body in replacement main")
    replacements = {}
    for path in inputs:
        candidate, bodies, _ = read_leb1(path)
        if (candidate.jd_start, candidate.jd_end) != (header.jd_start, header.jd_end):
            raise ValueError("replacement files must have the same global date range as the main")
        if candidate.body_count != len(bodies):
            raise ValueError("duplicate body in replacement partial")
        for body in bodies:
            identifier = body.entry.body_id
            if identifier in replacements:
                raise ValueError(f"duplicate replacement body: {identifier}")
            replacements[identifier] = body
    if replacements.keys() != expected_bodies:
        missing = sorted(expected_bodies - replacements.keys())
        unexpected = sorted(replacements.keys() - expected_bodies)
        raise ValueError(f"replacement inventory differs: missing={missing}, unexpected={unexpected}")
    retained = {body.entry.body_id: body for body in original}
    retained.update(replacements)
    bodies = [retained[identifier] for identifier in sorted(retained)]

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output.name}-replace-", dir=output.parent) as name:
        staged = Path(name) / "verified.leb"
        write_leb1(staged, bodies, header, auxiliary)
        if permissions is not None:
            staged.chmod(permissions)
        _check_replacement_paths(main, inputs, output, overwrite)
        check_input_hashes(identities)
        if overwrite:
            os.replace(staged, output)
        else:
            os.link(staged, output)


def _section_digest(section: SectionData) -> bytes:
    """Hash one bounded section without allocating the complete payload."""
    digest = hashlib.sha256()
    with section.path.open("rb") as source:
        source.seek(section.offset)
        remaining = section.size
        while remaining:
            block = source.read(min(remaining, 1024 * 1024))
            if not block:
                raise ValueError("truncated auxiliary source")
            digest.update(block)
            remaining -= len(block)
    return digest.digest()
