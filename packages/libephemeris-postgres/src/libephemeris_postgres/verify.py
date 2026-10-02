# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Streaming verification of imported LEB2 artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from libephemeris.leb_reader import open_leb

from .config import admin_dsn
from .importer import _artifact_key, _digest, _metadata


def verify_files(
    paths: list[str | Path], *, dataset_id: str, dsn: str | None = None
) -> None:
    """Verify checksums, metadata, and coefficient pages byte for byte."""

    try:
        import psycopg

        conn = psycopg.connect(admin_dsn(dsn), autocommit=False)
    except Exception:
        from libephemeris import CoefficientSourceError

        raise CoefficientSourceError("PostgreSQL provisioning failed") from None
    try:
        with conn.cursor() as cur:
            for raw_path in paths:
                path = Path(raw_path)
                _name, _group = _artifact_key(path, _tier_for_dataset(cur, dataset_id))
                digest = _digest(path)
                cur.execute(
                    "SELECT artifact_no, sha256 FROM libephemeris.artifacts "
                    "WHERE dataset_id = %s AND name = %s",
                    (dataset_id, path.name),
                )
                row = cur.fetchone()
                if row is None or row[1] != digest:
                    raise ValueError(f"checksum mismatch for {path.name}")
                _verify_one(cur, dataset_id, int(row[0]), path)
    finally:
        conn.close()


def _tier_for_dataset(cur: Any, dataset_id: str) -> str:
    cur.execute(
        "SELECT tier FROM libephemeris.datasets WHERE dataset_id = %s", (dataset_id,)
    )
    row = cur.fetchone()
    if row is None:
        raise ValueError("dataset not found")
    return str(row[0])


def _verify_one(cur: Any, dataset: str, artifact_no: int, path: Path) -> None:
    reader = open_leb(str(path))
    try:
        cur.execute(
            "SELECT jd_start, jd_end, delta_t_jd, delta_t_val "
            "FROM libephemeris.artifacts WHERE dataset_id=%s AND artifact_no=%s",
            (dataset, artifact_no),
        )
        artifact = cur.fetchone()
        if (
            artifact is None
            or float(artifact[0]).hex() != float(reader.header_jd_range[0]).hex()
            or float(artifact[1]).hex() != float(reader.header_jd_range[1]).hex()
            or [float(v).hex() for v in artifact[2]]
            != [float(v).hex() for v in reader.delta_t_table[0]]
            or [float(v).hex() for v in artifact[3]]
            != [float(v).hex() for v in reader.delta_t_table[1]]
        ):
            raise ValueError(f"metadata mismatch for {path.name}")
        entries = list(reader.bodies.items())
        if (
            reader.nutation_header is not None
            and reader.nutation_header.segment_count > 0
        ):
            entries.append((-1, reader.nutation_header))
        cur.execute(
            "SELECT body_id FROM libephemeris.series WHERE dataset_id=%s AND artifact_no=%s",
            (dataset, artifact_no),
        )
        if {row[0] for row in cur.fetchall()} != {body_id for body_id, _ in entries}:
            raise ValueError(f"series inventory mismatch for {path.name}")
        for body_id, entry in entries:
            cur.execute(
                "SELECT coord_type, segment_count, jd_start, jd_end, interval_days, degree, components "
                "FROM libephemeris.series WHERE dataset_id=%s AND artifact_no=%s AND body_id=%s",
                (dataset, artifact_no, body_id),
            )
            row = cur.fetchone()
            expected = _metadata(reader, body_id, entry)
            if (
                row is None
                or tuple(row[:2]) != expected[:2]
                or any(
                    float(db).hex() != float(want).hex()
                    for db, want in zip(row[2:5], expected[2:5])
                )
                or tuple(row[5:]) != expected[5:]
            ):
                raise ValueError(f"series metadata mismatch for {path.name}/{body_id}")
            with cur.connection.cursor(name="verify_pages", binary=True) as pages:
                pages.execute(
                    "SELECT page_no, coeffs FROM libephemeris.pages WHERE dataset_id=%s AND artifact_no=%s "
                    "AND body_id=%s ORDER BY page_no",
                    (dataset, artifact_no, body_id),
                )
                for page_no, _first, payload in reader.iter_segment_pages(
                    body_id, page_size=32
                ):
                    stored = pages.fetchone()
                    if (
                        stored is None
                        or stored[0] != page_no
                        or bytes(stored[1]) != payload
                    ):
                        raise ValueError(
                            f"page mismatch for {path.name}/{body_id}/{page_no}"
                        )
                if pages.fetchone() is not None:
                    raise ValueError(
                        f"extra coefficient pages for {path.name}/{body_id}"
                    )
        cur.execute(
            "SELECT star_id, ra_j2000, dec_j2000, pm_ra, pm_dec, parallax, rv, magnitude "
            "FROM libephemeris.stars WHERE dataset_id=%s AND artifact_no=%s ORDER BY star_id",
            (dataset, artifact_no),
        )
        expected_stars = reader.stars
        rows = cur.fetchall()
        if {row[0] for row in rows} != set(expected_stars):
            raise ValueError(f"star inventory mismatch for {path.name}")
        for row in rows:
            star = expected_stars[row[0]]
            expected = (
                star.ra_j2000,
                star.dec_j2000,
                star.pm_ra,
                star.pm_dec,
                star.parallax,
                star.rv,
                star.magnitude,
            )
            if any(float(a).hex() != float(b).hex() for a, b in zip(row[1:], expected)):
                raise ValueError(f"star metadata mismatch for {path.name}")
    finally:
        reader.close()
