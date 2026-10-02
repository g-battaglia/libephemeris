# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Opt-in integration checks for a disposable PostgreSQL instance."""

from __future__ import annotations

import os
import uuid
from pathlib import Path

from libephemeris.leb_reader import open_leb

import pytest

from libephemeris_postgres import open_tier
from libephemeris_postgres.importer import import_files
from libephemeris_postgres.verify import verify_files
from libephemeris_postgres.__main__ import _schema

psycopg = pytest.importorskip("psycopg")


pytestmark = pytest.mark.pg_integration

_DSN = os.environ.get("LIBEPHEMERIS_TEST_PG_URL")
if not _DSN:
    pytestmark = [
        pytest.mark.pg_integration,
        pytest.mark.skip(reason="LIBEPHEMERIS_TEST_PG_URL is not set"),
    ]

_ROOT = Path(__file__).resolve().parents[3]
_CORE = _ROOT / "libephemeris" / "data" / "leb2" / "base_core.leb2"


def _four_artifacts(tmp_path: Path) -> list[Path]:
    """Use the bundled core bytes under each required group name."""

    paths = []
    for group in ("core", "asteroids", "exotics", "apogee"):
        target = tmp_path / f"base_{group}.leb2"
        target.write_bytes(_CORE.read_bytes())
        paths.append(target)
    return paths


@pytest.mark.skipif(not _DSN, reason="LIBEPHEMERIS_TEST_PG_URL is not set")
def test_import_resume_verify_and_corruption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exercise the complete artifact lifecycle against an explicit test DSN."""

    assert _CORE.exists()
    monkeypatch.setenv("LIBEPHEMERIS_PG_ADMIN_URL", _DSN)
    _schema(_DSN)
    paths = _four_artifacts(tmp_path)
    dataset = str(uuid.uuid4())
    import_files(paths, tier="base", dsn=_DSN, dataset_id=dataset)
    with psycopg.connect(_DSN, autocommit=True) as connection:
        with connection.cursor() as cursor:
            cursor.execute(
                "SELECT count(*) FROM libephemeris.artifacts WHERE dataset_id = %s",
                (dataset,),
            )
            artifact_count = cursor.fetchone()[0]
            cursor.execute(
                "SELECT count(*) FROM libephemeris.pages WHERE dataset_id = %s",
                (dataset,),
            )
            page_count = cursor.fetchone()[0]
    import_files(paths, tier="base", dsn=_DSN, resume=dataset)
    with psycopg.connect(_DSN, autocommit=True) as connection:
        with connection.cursor() as cursor:
            cursor.execute(
                "SELECT count(*) FROM libephemeris.artifacts WHERE dataset_id = %s",
                (dataset,),
            )
            assert cursor.fetchone()[0] == artifact_count
            cursor.execute(
                "SELECT count(*) FROM libephemeris.pages WHERE dataset_id = %s",
                (dataset,),
            )
            assert cursor.fetchone()[0] == page_count
    verify_files(paths, dataset_id=dataset, dsn=_DSN)
    monkeypatch.setenv("LIBEPHEMERIS_PG_URL", _DSN)
    monkeypatch.setenv("LIBEPHEMERIS_PG_DATASET_BASE", dataset)
    sources = open_tier("base")
    assert len(sources) == 4

    file_reader = open_leb(str(_CORE))
    try:
        source = sources[0]
        for body_id, entry in file_reader.bodies.items():
            for idx in (0, entry.segment_count // 2, entry.segment_count - 1):
                assert source.eval_body(
                    body_id, entry.jd_start + idx * entry.interval_days
                ) == file_reader.eval_body(
                    body_id, entry.jd_start + idx * entry.interval_days
                )
        for jd, expected in zip(*file_reader.delta_t_table):
            assert source.delta_t(jd) == expected
        if file_reader.nutation_header is not None:
            for jd in (
                file_reader.nutation_header.jd_start,
                file_reader.nutation_header.jd_end,
            ):
                assert source.eval_nutation(jd) == file_reader.eval_nutation(jd)
    finally:
        file_reader.close()

    with psycopg.connect(_DSN, autocommit=True) as connection:
        with connection.cursor() as cursor:
            cursor.execute(
                "UPDATE libephemeris.pages SET coeffs = set_byte(coeffs, 0, "
                "get_byte(coeffs, 0) # 255) WHERE dataset_id = %s AND body_id = 0 "
                "AND page_no = 0",
                (dataset,),
            )
    with pytest.raises(ValueError, match="page mismatch"):
        verify_files(paths, dataset_id=dataset, dsn=_DSN)


def test_runtime_role_is_read_only() -> None:
    """Document the optional runtime-role check without requiring credentials."""

    runtime_dsn = os.environ.get("LIBEPHEMERIS_TEST_PG_RUNTIME_URL")
    if not runtime_dsn:
        pytest.skip("LIBEPHEMERIS_TEST_PG_RUNTIME_URL is not set")
    with psycopg.connect(runtime_dsn, autocommit=True) as connection:
        with connection.cursor() as cursor:
            with pytest.raises(psycopg.errors.ReadOnlySqlTransaction):
                cursor.execute(
                    "CREATE TABLE libephemeris._provider_write_probe (id integer)"
                )
            cursor.execute("SET default_transaction_read_only=off")
            with pytest.raises(psycopg.errors.InsufficientPrivilege):
                cursor.execute("INSERT INTO libephemeris.schema_version VALUES (1)")
