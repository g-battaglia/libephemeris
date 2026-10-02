# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Explicit offline import of native LEB artifacts into PostgreSQL.

Provenance:
    Lossless storage conversion of project-native artifacts with their registered
    JPL/IAU provenance. Existing Zstandard-based LEB2 decompression is reused;
    IEEE-754 coefficients are never resampled or refitted. The storage contract
    and publication policy are documented in docs/db-backend.md.
"""

from __future__ import annotations

import hashlib
import json
from contextlib import ExitStack
from importlib.resources import files
from pathlib import Path
from typing import Any, Iterator
from uuid import UUID

from ..exceptions import DBDataError, DBError
from ..leb_compression import decompress_body
from ..leb_format import (
    HEADER_SIZE,
    LEB2_MAGIC,
    SECTION_CHEBYSHEV,
    SECTION_COMPRESSED_CHEBYSHEV,
)
from ..leb_reader import open_leb
from .contract import (
    NUTATION_ID,
    Series,
    decode_segment,
    segment_byte_count,
    series_fields,
    validate_series,
)


class Artifact:
    """Own an offline file reader and expose normalized records for import."""

    def __init__(self, path: Path, reader: Any) -> None:
        """Validate section bounds and derive source identity.

        Args:
            path: Source artifact, opened only by the offline importer.
            reader: Existing LEB1 or LEB2 parser.

        Raises:
            DBDataError: Section directory refers outside the artifact.
        """
        self.path = path
        self.reader = reader
        if reader._header.body_count != len(reader._bodies):
            raise DBDataError("Duplicate or missing body entries in source artifact")
        if reader._header.section_count != len(reader._sections):
            raise DBDataError("Duplicate section identifiers in source artifact")
        for section in reader._sections.values():
            if section.offset < 0 or section.offset + section.size > len(reader._mm):
                raise DBDataError("LEB section outside artifact bounds")
        with path.open("rb") as source:
            self.sha256 = hashlib.file_digest(source, "sha256").hexdigest()

    def manifest_record(self) -> dict:
        """Describe the artifact without exposing filesystem-specific paths.

        Returns:
            JSON-compatible source hash, header, and section-directory fields.
        """
        source_header = self.reader._header
        header = {
            "magic": source_header.magic.decode("ascii"),
            "version": source_header.version,
            "section_count": source_header.section_count,
            "body_count": source_header.body_count,
            "jd_start": source_header.jd_start,
            "jd_end": source_header.jd_end,
            "generation_epoch": source_header.generation_epoch,
            "flags": source_header.flags,
        }
        sections = []
        for section in self.reader._sections.values():
            sections.append(
                {
                    "section_id": section.section_id,
                    "offset": section.offset,
                    "size": section.size,
                }
            )
        return {
            "name": self.path.name,
            "sha256": self.sha256,
            "header": header,
            "header_hex": bytes(self.reader._mm[:HEADER_SIZE]).hex(),
            "sections": sections,
        }

    def body_series(self) -> Iterator[Series]:
        """Normalize source body headers into portable coefficient metadata.

        Yields:
            One validated Series for each source body.
        """
        for body in self.reader._bodies.values():
            series = Series(
                body.body_id,
                body.coord_type,
                body.segment_count,
                body.jd_start,
                body.jd_end,
                body.interval_days,
                body.degree,
                body.components,
            )
            validate_series(series)
            if series.body_id == NUTATION_ID:
                raise DBDataError("Source body identifier conflicts with nutation")
            yield series

    def body_segments(self, series: Series) -> Iterator[tuple[int, bytes]]:
        """Decode coefficients in bounded chunks rather than retaining a dataset.

        Args:
            series: Validated body metadata belonging to this artifact.

        Yields:
            Zero-based global segment index and exact little-endian bytes.

        Raises:
            DBDataError: Payload or chunk indices do not match body metadata.
        """
        reader = self.reader
        body = reader._bodies[series.body_id]
        size = segment_byte_count(series)
        if reader._header.magic == LEB2_MAGIC:
            if reader._chunked:
                next_index = 0
                for chunk in reader._chunk_index[series.body_id]:
                    if chunk.segment_start != next_index:
                        raise DBDataError(
                            "Source chunks do not cover consecutive segments"
                        )
                    expected_size = chunk.segment_count * size
                    if chunk.uncompressed_size != expected_size:
                        raise DBDataError("Invalid source chunk length")
                    compressed = self._slice(chunk.blob_offset, chunk.compressed_size)
                    decoded = decompress_body(
                        compressed,
                        expected_size,
                        chunk.segment_count,
                        series.degree,
                        series.components,
                    )
                    for local_index in range(chunk.segment_count):
                        offset = local_index * size
                        yield next_index + local_index, decoded[offset : offset + size]
                    next_index += chunk.segment_count
                if next_index != series.segment_count:
                    raise DBDataError("Source chunk inventory is incomplete")
            else:
                expected_size = series.segment_count * size
                if body.uncompressed_size != expected_size:
                    raise DBDataError("Invalid monolithic source length")
                compressed = self._slice(body.data_offset, body.compressed_size)
                decoded = decompress_body(
                    compressed,
                    expected_size,
                    series.segment_count,
                    series.degree,
                    series.components,
                )
                for index in range(series.segment_count):
                    yield index, decoded[index * size : (index + 1) * size]
        else:
            for index in range(series.segment_count):
                yield index, self._slice(body.data_offset + index * size, size)

    def _slice(self, offset: int, size: int) -> bytes:
        """Read an exact source range and reject silent mmap truncation.

        Args:
            offset: Absolute source byte offset.
            size: Required number of bytes.

        Returns:
            Detached source bytes.

        Raises:
            DBDataError: Range is invalid or truncated.
        """
        if offset < 0 or size < 0 or offset + size > len(self.reader._mm):
            raise DBDataError("Truncated LEB coefficient or section payload")
        return bytes(self.reader._mm[offset : offset + size])

    def auxiliary_sections(self) -> Iterator[tuple[int, bytes]]:
        """Preserve non-coefficient sections, including reserved opaque data.

        Yields:
            Section identifier and original payload. Coefficients are stored
            separately as decoded segments, never duplicated as large blobs.
        """
        for section in self.reader._sections.values():
            if section.section_id not in (
                SECTION_CHEBYSHEV,
                SECTION_COMPRESSED_CHEBYSHEV,
            ):
                yield section.section_id, self._slice(section.offset, section.size)


def provision_schema(connection: Any) -> None:
    """Install schema v1 explicitly using a migration-capable connection.

    Args:
        connection: PostgreSQL connection with schema-creation privileges.
    """
    sql = files("libephemeris.db").joinpath("schema_v1.sql").read_text()
    connection.execute(sql)


def _copy_records(connection: Any, table: str, columns: str, records: Iterator) -> None:
    """Stream trusted importer records through PostgreSQL COPY.

    Args:
        connection: Provisioning connection in the import transaction.
        table: Internal constant table name, never user-controlled.
        columns: Internal constant list of columns.
        records: Iterator producing native record tuples.
    """
    with connection.cursor().copy(
        f"COPY libephemeris.{table} ({columns}) FROM STDIN"
    ) as copy:
        for record in records:
            copy.write_row(record)


def _segment_records(
    dataset_id: str,
    series: Series,
    segments: Iterator[tuple[int, bytes]],
) -> Iterator[tuple[str, int, int, bytes]]:
    """Validate a consecutive segment stream before bulk insertion.

    Args:
        dataset_id: Immutable dataset being constructed.
        series: Metadata describing the source stream.
        segments: Index/payload pairs decoded without approximation.

    Yields:
        SQL records containing dataset, body, index and coefficient bytes.

    Raises:
        DBDataError: Stream is incomplete, unordered or contains invalid values.
    """
    expected_index = 0
    for index, payload in segments:
        if index != expected_index:
            raise DBDataError("Source segments are not consecutive")
        decode_segment(series, payload)
        yield dataset_id, series.body_id, index, payload
        expected_index += 1
    if expected_index != series.segment_count:
        raise DBDataError("Source segment inventory is incomplete")


def _import_body(
    connection: Any, dataset_id: str, artifact: Artifact, series: Series
) -> None:
    """Insert one header and its complete coefficient inventory.

    Args:
        connection: Provisioning connection in the import transaction.
        dataset_id: Immutable version UUID being constructed.
        artifact: Offline source owner.
        series: Validated metadata to import.
    """
    connection.execute(
        "INSERT INTO libephemeris.series VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)",
        (dataset_id, *series_fields(series)),
    )
    records = _segment_records(dataset_id, series, artifact.body_segments(series))
    _copy_records(
        connection, "segments", "dataset_id,body_id,segment_index,payload", records
    )


def _import_auxiliary(
    connection: Any, dataset_id: str, artifacts: list[Artifact]
) -> None:
    """Normalize auxiliary records and preserve every non-coefficient section.

    Companion files may repeat identical auxiliary sections. Conflicting
    realizations are rejected instead of silently choosing a different model.

    Args:
        connection: Provisioning connection in the import transaction.
        dataset_id: Version being constructed.
        artifacts: Ordered source owners from one tier.

    Raises:
        DBDataError: Repeated auxiliary data disagrees across companions.
    """
    nutation_identity = None
    delta_t_identity = None
    star_identity = None
    for artifact_index, artifact in enumerate(artifacts):
        reader = artifact.reader
        for section_id, payload in artifact.auxiliary_sections():
            connection.execute(
                "INSERT INTO libephemeris.sections VALUES (%s,%s,%s,%s)",
                (dataset_id, artifact_index, section_id, payload),
            )
        nutation = reader._nutation
        if nutation is not None and nutation.segment_count > 0:
            series = Series(
                NUTATION_ID,
                0,
                nutation.segment_count,
                nutation.jd_start,
                nutation.jd_end,
                nutation.interval_days,
                nutation.degree,
                nutation.components,
            )
            validate_series(series)
            segment_size = segment_byte_count(series)
            payload = artifact._slice(
                reader._nutation_data_offset,
                series.segment_count * segment_size,
            )
            identity = (series, hashlib.sha256(payload).digest())
            if nutation_identity is not None and identity != nutation_identity:
                raise DBDataError("Companion nutation sections disagree")
            if nutation_identity is None:
                connection.execute(
                    "INSERT INTO libephemeris.series VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)",
                    (dataset_id, *series_fields(series)),
                )
                segments = (
                    (
                        index,
                        payload[index * segment_size : (index + 1) * segment_size],
                    )
                    for index in range(series.segment_count)
                )
                records = _segment_records(dataset_id, series, segments)
                _copy_records(
                    connection,
                    "segments",
                    "dataset_id,body_id,segment_index,payload",
                    records,
                )
                nutation_identity = identity
        if reader._delta_t_jds:
            points = list(zip(reader._delta_t_jds, reader._delta_t_vals))
            if delta_t_identity is not None and points != delta_t_identity:
                raise DBDataError("Companion Delta-T sections disagree")
            if delta_t_identity is None:
                _copy_records(
                    connection,
                    "delta_t",
                    "dataset_id,jd,days",
                    ((dataset_id, jd, days) for jd, days in points),
                )
                delta_t_identity = points
        if reader._stars:
            entries = []
            for star in reader._stars.values():
                entries.append(
                    (
                        star.star_id,
                        star.ra_j2000,
                        star.dec_j2000,
                        star.pm_ra,
                        star.pm_dec,
                        star.parallax,
                        star.rv,
                        star.magnitude,
                    )
                )
            if star_identity is not None and entries != star_identity:
                raise DBDataError("Companion star sections disagree")
            if star_identity is None:
                _copy_records(
                    connection,
                    "stars",
                    "dataset_id,star_id,ra,dec,pm_ra,pm_dec,parallax,rv,magnitude",
                    ((dataset_id, *entry) for entry in entries),
                )
                star_identity = entries


def import_artifacts(
    connection: Any,
    dataset_id: str,
    paths: list[Path],
    tier: str,
) -> None:
    """Publish a complete immutable dataset in one offline transaction.

    Args:
        connection: Separate read-write provisioning connection.
        dataset_id: Explicit UUID for the new version.
        paths: Explicit LEB1/LEB2 artifacts from a single tier.
        tier: Dataset tier: base, medium or extended.

    Raises:
        DBDataError: Artifacts conflict or an existing identity differs.
        ValueError: UUID or tier is invalid.
    """
    dataset_id = str(UUID(dataset_id))
    if tier not in ("base", "medium", "extended") or not paths:
        raise ValueError("Import requires artifacts and a valid explicit tier")
    with ExitStack() as stack:
        artifacts = []
        for path in paths:
            encoded_tiers = set(path.stem.split("_")) & {"base", "medium", "extended"}
            if encoded_tiers and encoded_tiers != {tier}:
                raise DBDataError("Cannot mix artifact tiers in one DB dataset")
            try:
                reader = open_leb(str(path))
            except OSError as error:
                raise DBDataError(
                    "Coefficient source artifact is unavailable or unreadable"
                ) from error
            stack.callback(reader.close)
            artifacts.append(Artifact(path, reader))
        try:
            manifest = {
                "artifacts": [artifact.manifest_record() for artifact in artifacts]
            }
        except OSError as error:
            raise DBDataError(
                "Coefficient source artifact is unavailable or unreadable"
            ) from error
        with connection.transaction():
            # Serializes two provisioners targeting the same immutable UUID.
            connection.execute(
                "SELECT pg_advisory_xact_lock(hashtextextended(%s, 0))", (dataset_id,)
            )
            existing = connection.execute(
                "SELECT tier, published, manifest FROM libephemeris.datasets WHERE dataset_id=%s",
                (dataset_id,),
            ).fetchone()
            if existing is not None:
                if existing == (tier, True, manifest):
                    return
                raise DBDataError(
                    "Dataset UUID already belongs to different/incomplete artifacts"
                )
            connection.execute(
                "INSERT INTO libephemeris.datasets (dataset_id,tier,manifest) VALUES (%s,%s,%s::jsonb)",
                (dataset_id, tier, json.dumps(manifest)),
            )
            seen_bodies = set()
            for artifact in artifacts:
                for series in artifact.body_series():
                    if series.body_id in seen_bodies:
                        raise DBDataError(
                            "Duplicate body channel across source artifacts"
                        )
                    _import_body(connection, dataset_id, artifact, series)
                    seen_bodies.add(series.body_id)
            _import_auxiliary(connection, dataset_id, artifacts)
            if not seen_bodies:
                raise DBDataError("Cannot publish a dataset without body channels")
            connection.execute(
                "UPDATE libephemeris.datasets SET published=true WHERE dataset_id=%s",
                (dataset_id,),
            )
        # Bulk COPY can outpace autovacuum's statistics refresh. Without fresh
        # estimates, a new large dataset can get a sequential-scan join plan
        # and hit the sealed runtime statement limit. This is an explicit
        # provisioning action, never a calculation-time write or migration.
        connection.execute(
            "ANALYZE libephemeris.datasets, libephemeris.series, libephemeris.segments"
        )


def connect_provisioner(dsn: str) -> Any:
    """Open an explicit writer connection without exposing credentials.

    Args:
        dsn: Provisioning PostgreSQL connection string.

    Returns:
        A writer connection, to be used as a context manager by the caller.

    Raises:
        DBError: Optional driver is missing or connection fails.
    """
    from .policy import require_database

    require_database()
    try:
        import psycopg
    except ImportError:
        raise DBError("Install libephemeris[postgres] to import datasets") from None
    try:
        return psycopg.connect(dsn, autocommit=True, connect_timeout=10)
    except psycopg.Error:
        raise DBError("Could not connect to the provisioning database") from None
