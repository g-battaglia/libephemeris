# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disposable PostgreSQL checks for statement guards and publication locks."""

from __future__ import annotations

import os
import struct
from uuid import uuid4

import pytest

import libephemeris as ephe
from libephemeris.db.importer import (
    connect_provisioner,
    import_artifacts,
    provision_schema,
)
from .artifacts import make_artifacts


@pytest.fixture
def database(tmp_path):
    dsn = os.environ.get("LIBEPHEMERIS_TEST_DB_URL")
    if not dsn:
        pytest.skip("Requires an explicitly selected disposable PostgreSQL database")
    pytest.importorskip("psycopg")
    policy = ephe.get_configured_network_policy()
    ephe.set_network_policy("allow")
    try:
        with connect_provisioner(dsn) as connection:
            with connection.transaction():
                provision_schema(connection)
                provision_schema(connection)
            published = str(uuid4())
            import_artifacts(
                connection, published, [make_artifacts(tmp_path)[0]], "base"
            )
            yield connection, dsn, published
    finally:
        ephe.set_network_policy(policy)


@pytest.mark.db_integration
@pytest.mark.parametrize(
    "table",
    [
        "schema_version",
        "datasets",
        "series",
        "segments",
        "delta_t",
        "stars",
        "sections",
    ],
)
def test_truncate_cannot_remove_scientific_publications(database, table):
    import psycopg

    connection, _, published = database
    with pytest.raises(psycopg.Error, match="cannot be truncated"):
        connection.execute(f"TRUNCATE libephemeris.{table} CASCADE")
    assert connection.execute(
        "SELECT published FROM libephemeris.datasets WHERE dataset_id=%s", (published,)
    ).fetchone() == (True,)


def _unpublished(connection, segments=2):
    dataset = str(uuid4())
    connection.execute(
        "INSERT INTO libephemeris.datasets VALUES (%s,'base',false,'{}')", (dataset,)
    )
    connection.execute(
        "INSERT INTO libephemeris.series VALUES (%s,0,0,%s,0,%s,1,0,3)",
        (dataset, segments, float(segments)),
    )
    return dataset


@pytest.mark.db_integration
def test_segment_copy_checks_publication_once_per_statement_and_rolls_back(database):
    import psycopg

    connection, _, published = database
    triggers = dict(
        connection.execute(
            "SELECT tgname,tgtype FROM pg_trigger WHERE tgrelid='libephemeris.segments'::regclass AND NOT tgisinternal"
        ).fetchall()
    )
    assert triggers["immutable_insert_segments"] == 4  # AFTER INSERT, statement
    assert triggers["immutable_segments"] == 27  # BEFORE UPDATE/DELETE, row
    unpublished = _unpublished(connection)
    with pytest.raises(psycopg.Error, match="immutable"):
        with connection.transaction(), connection.cursor() as cursor:
            with cursor.copy("COPY libephemeris.segments FROM STDIN") as copy:
                copy.write_row((unpublished, 0, 0, struct.pack("<3d", 1.0, 2.0, 3.0)))
                copy.write_row((published, 0, 2, struct.pack("<3d", 1.0, 2.0, 3.0)))
    assert connection.execute(
        "SELECT count(*) FROM libephemeris.segments WHERE dataset_id=%s", (unpublished,)
    ).fetchone() == (0,)


@pytest.mark.db_integration
def test_bulk_insert_locks_parent_until_publication_is_safe(database):
    import psycopg

    connection, dsn, _ = database
    count = 1000
    unpublished = _unpublished(connection, count)
    payload = struct.pack("<3d", 1.0, 2.0, 3.0)
    with connect_provisioner(dsn) as publisher:
        with connection.transaction(), connection.cursor() as cursor:
            with cursor.copy("COPY libephemeris.segments FROM STDIN") as copy:
                for index in range(count):
                    copy.write_row((unpublished, 0, index, payload))
            with pytest.raises(psycopg.errors.LockNotAvailable):
                with publisher.transaction():
                    publisher.execute("SET LOCAL lock_timeout='100ms'")
                    publisher.execute(
                        "UPDATE libephemeris.datasets SET published=true WHERE dataset_id=%s",
                        (unpublished,),
                    )
        publisher.execute(
            "UPDATE libephemeris.datasets SET published=true WHERE dataset_id=%s",
            (unpublished,),
        )
    assert connection.execute(
        "SELECT count(*) FROM libephemeris.segments WHERE dataset_id=%s", (unpublished,)
    ).fetchone() == (count,)
    with pytest.raises(psycopg.Error, match="immutable"):
        connection.execute(
            "UPDATE libephemeris.segments SET payload=%s WHERE dataset_id=%s AND segment_index=0",
            (payload, unpublished),
        )
