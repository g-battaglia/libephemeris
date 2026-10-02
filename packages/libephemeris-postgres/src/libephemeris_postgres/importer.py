# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Resumable, artifact-at-a-time LEB2 importer."""

from __future__ import annotations

import hashlib
import math
import re
import uuid
from pathlib import Path
from typing import Any

from libephemeris.leb_groups import LEB2_GROUPS
from libephemeris.leb_reader import open_leb

from .config import admin_dsn

_ARTIFACT_RE = re.compile(r"^(base|medium|extended)_([a-z]+)\.leb2$")


def _own_error(message: str) -> Exception:
    from libephemeris import CoefficientSourceError

    return CoefficientSourceError(message)


def _artifact_key(path: Path, tier: str) -> tuple[str, str]:
    match = _ARTIFACT_RE.fullmatch(path.name)
    if match is None or match.group(1) != tier or match.group(2) not in LEB2_GROUPS:
        raise _own_error(f"invalid artifact name: {path.name}")
    return path.name, match.group(2)


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _metadata(reader: Any, body_id: int, entry: Any) -> tuple[Any, ...]:
    values = (
        int(getattr(entry, "coord_type", 0)),
        int(getattr(entry, "segment_count")),
        float(getattr(entry, "jd_start")),
        float(getattr(entry, "jd_end")),
        float(getattr(entry, "interval_days")),
        int(getattr(entry, "degree")),
        int(getattr(entry, "components")),
    )
    if values[1] <= 0 or values[3] <= values[2] or values[4] <= 0:
        raise _own_error(f"invalid metadata for body {body_id}")
    if not all(math.isfinite(value) for value in values[2:5]):
        raise _own_error(f"non-finite metadata for body {body_id}")
    return values


def import_files(
    paths: list[str | Path],
    *,
    tier: str,
    dsn: str | None = None,
    dataset_id: str | None = None,
    resume: str | None = None,
    libephemeris_version: str = "unknown",
) -> str:
    """Import artifacts, committing each artifact in its own transaction."""

    if tier not in {"base", "medium", "extended"}:
        raise _own_error(f"invalid tier: {tier}")
    dataset = resume or dataset_id or str(uuid.uuid4())
    try:
        import psycopg

        conn = psycopg.connect(admin_dsn(dsn), autocommit=False)
    except Exception:
        raise _own_error("PostgreSQL provisioning failed") from None
    try:
        with conn.cursor() as cur:
            # Session lock spans artifact commits and prevents concurrent resume.
            cur.execute("SELECT pg_advisory_lock(hashtextextended(%s, 0))", (dataset,))
            if resume:
                cur.execute(
                    "SELECT tier, complete FROM libephemeris.datasets WHERE dataset_id=%s",
                    (dataset,),
                )
                existing_dataset = cur.fetchone()
                if existing_dataset is None or existing_dataset[0] != tier:
                    raise _own_error(
                        "Resume dataset is absent or belongs to another tier"
                    )
                if existing_dataset[1]:
                    return dataset
            if not resume:
                cur.execute(
                    "INSERT INTO libephemeris.datasets "
                    "(dataset_id, tier, libephemeris_version) VALUES (%s, %s, %s)",
                    (dataset, tier, libephemeris_version),
                )
            ordered = sorted(
                paths,
                key=lambda value: (
                    not Path(value).name.endswith("_core.leb2"),
                    Path(value).name,
                ),
            )
            for raw_path in ordered:
                path = Path(raw_path)
                name, group = _artifact_key(path, tier)
                artifact_no = LEB2_GROUPS.index(group)
                sha = _digest(path)
                cur.execute(
                    "SELECT artifact_no, sha256 FROM libephemeris.artifacts "
                    "WHERE dataset_id = %s AND name = %s",
                    (dataset, name),
                )
                existing = cur.fetchone()
                if existing:
                    if existing[1] != sha:
                        raise _own_error(f"artifact {name} has a different checksum")
                    continue
                reader = open_leb(str(path))
                try:
                    _import_one(cur, dataset, artifact_no, name, sha, reader)
                finally:
                    reader.close()
                conn.commit()
            cur.execute(
                "SELECT count(*) FROM libephemeris.artifacts WHERE dataset_id = %s",
                (dataset,),
            )
            count_row = cur.fetchone()
            if count_row is not None and count_row[0] == len(LEB2_GROUPS):
                cur.execute(
                    "UPDATE libephemeris.datasets SET complete = true WHERE dataset_id = %s",
                    (dataset,),
                )
                cur.execute(
                    "ANALYZE libephemeris.datasets, libephemeris.artifacts, "
                    "libephemeris.series, libephemeris.pages"
                )
            conn.commit()
    except psycopg.Error:
        conn.rollback()
        raise _own_error("PostgreSQL provisioning failed") from None
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()
    return dataset


def _copy_pages(cur: Any, rows: Any) -> None:
    """Use COPY when available, while keeping the writer easy to test."""

    with cur.copy(
        "COPY libephemeris.pages "
        "(dataset_id, artifact_no, body_id, page_no, coeffs) FROM STDIN (FORMAT BINARY)"
    ) as copy:
        copy.set_types(["uuid", "int2", "int4", "int4", "bytea"])
        for dataset, artifact, body, page, payload in rows:
            copy.write_row((uuid.UUID(dataset), artifact, body, page, payload))


def _import_one(
    cur: Any, dataset: str, artifact_no: int, name: str, sha: str, reader: Any
) -> None:
    jd_start, jd_end = reader.header_jd_range
    delta_jd, delta_values = reader.delta_t_table
    cur.execute(
        "INSERT INTO libephemeris.artifacts "
        "(dataset_id, artifact_no, name, sha256, jd_start, jd_end, delta_t_jd, delta_t_val) "
        "VALUES (%s,%s,%s,%s,%s,%s,%s,%s)",
        (
            dataset,
            artifact_no,
            name,
            sha,
            jd_start,
            jd_end,
            list(delta_jd),
            list(delta_values),
        ),
    )
    entries = list(reader.bodies.items())
    if reader.nutation_header is not None and reader.nutation_header.segment_count > 0:
        entries.append((-1, reader.nutation_header))
    for body_id, entry in entries:
        values = _metadata(reader, body_id, entry)
        cur.execute(
            "INSERT INTO libephemeris.series "
            "(dataset_id, artifact_no, body_id, coord_type, segment_count, jd_start, jd_end, "
            "interval_days, degree, components, page_size) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            (dataset, artifact_no, body_id, *values, 32),
        )
        pages = (
            (dataset, artifact_no, body_id, page_no, payload)
            for page_no, _first, payload in reader.iter_segment_pages(
                body_id, page_size=32
            )
        )
        _copy_pages(cur, pages)
    for star_id, star in reader.stars.items():
        cur.execute(
            "INSERT INTO libephemeris.stars "
            "(dataset_id, artifact_no, star_id, ra_j2000, dec_j2000, pm_ra, pm_dec, "
            "parallax, rv, magnitude) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
            (
                dataset,
                artifact_no,
                star_id,
                star.ra_j2000,
                star.dec_j2000,
                star.pm_ra,
                star.pm_dec,
                star.parallax,
                star.rv,
                star.magnitude,
            ),
        )
