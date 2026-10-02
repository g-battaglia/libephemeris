# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Opt-in lifecycle and parity checks on disposable PostgreSQL storage."""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from libephemeris import CoefficientSourceError
from libephemeris.download import DATA_FILES
from libephemeris.leb2_reader import LEB2Reader
from libephemeris.leb_groups import LEB2_GROUPS

from libephemeris_postgres import open_tier
from libephemeris_postgres.__main__ import _schema
from libephemeris_postgres.importer import upload_files
from libephemeris_postgres.pool import reset_pool
from libephemeris_postgres.source import BLOCK_SIZE

psycopg = pytest.importorskip("psycopg")
_DSN = os.environ.get("LIBEPHEMERIS_TEST_PG_URL")
pytestmark = [
    pytest.mark.pg_integration,
    pytest.mark.skipif(not _DSN, reason="Disposable PG URL not set"),
]
_ROOT = Path(__file__).resolve().parents[3]


def test_upload_resume_verified_bytes_and_reader_parity(monkeypatch):
    paths = [_ROOT / "data" / "leb2" / f"medium_{group}.leb2" for group in LEB2_GROUPS]
    if not all(path.exists() for path in paths):
        pytest.skip("Real medium artifacts are not installed")
    _schema(_DSN)
    monkeypatch.setenv("LIBEPHEMERIS_PG_URL", _DSN)
    first = paths[0]
    sha = DATA_FILES[first.name]["sha256"]
    # Simulate an interrupted upload with its first block already committed.
    with psycopg.connect(_DSN) as conn:
        conn.execute(
            "INSERT INTO libephemeris.files (sha256,name,size) VALUES (%s,%s,%s)",
            (sha, first.name, first.stat().st_size),
        )
        with first.open("rb") as stream:
            conn.execute(
                "INSERT INTO libephemeris.blocks VALUES (%s,0,%s)",
                (sha, stream.read(BLOCK_SIZE)),
            )
    upload_files(paths, dsn=_DSN)
    upload_files(paths, dsn=_DSN)
    sources = open_tier("medium")
    try:
        for path in paths:
            remote = next(
                source for source in sources if source.path.endswith(path.name)
            )
            with LEB2Reader(str(path)) as local:
                for body, entry in local._bodies.items():
                    for jd in (
                        entry.jd_start,
                        (entry.jd_start + entry.jd_end) / 2,
                        entry.jd_end,
                    ):
                        assert remote.eval_body(body, jd) == local.eval_body(body, jd)
                if local.has_nutation():
                    for jd in (local._nutation.jd_start, local._nutation.jd_end):
                        assert remote.eval_nutation(jd) == local.eval_nutation(jd)
                assert remote._delta_t_jds == local._delta_t_jds
                assert remote._delta_t_vals == local._delta_t_vals
                assert remote._stars == local._stars
    finally:
        for source in sources:
            source.close()
        reset_pool()


def test_incomplete_pin_is_rejected(monkeypatch):
    sha = DATA_FILES["medium_core.leb2"]["sha256"]
    monkeypatch.setenv("LIBEPHEMERIS_PG_URL", _DSN)
    with psycopg.connect(_DSN) as conn:
        conn.execute(
            "UPDATE libephemeris.files SET complete=false WHERE sha256=%s", (sha,)
        )
    try:
        with pytest.raises(CoefficientSourceError):
            open_tier("medium")
    finally:
        with psycopg.connect(_DSN) as conn:
            conn.execute(
                "UPDATE libephemeris.files SET complete=true WHERE sha256=%s", (sha,)
            )
        reset_pool()


def test_corrupt_incomplete_upload_never_published(tmp_path, monkeypatch):
    import hashlib

    path = tmp_path / "base_core.leb2"
    data = b"test-byte-storage" * 5000
    path.write_bytes(data)
    sha = hashlib.sha256(data).hexdigest()
    monkeypatch.setitem(DATA_FILES, path.name, {"sha256": sha})
    with psycopg.connect(_DSN) as conn:
        conn.execute(
            "INSERT INTO libephemeris.files (sha256,name,size) VALUES (%s,%s,%s)",
            (sha, path.name, len(data)),
        )
        conn.execute(
            "INSERT INTO libephemeris.blocks VALUES (%s,0,%s)", (sha, b"x" * BLOCK_SIZE)
        )
    with pytest.raises(CoefficientSourceError, match="checksum"):
        upload_files([path], dsn=_DSN)
    with psycopg.connect(_DSN) as conn:
        assert conn.execute(
            "SELECT complete FROM libephemeris.files WHERE sha256=%s", (sha,)
        ).fetchone() == (False,)


def test_runtime_role_is_read_only():
    dsn = os.environ.get("LIBEPHEMERIS_TEST_PG_RUNTIME_URL")
    if not dsn:
        pytest.skip("Runtime-role URL not set")
    with psycopg.connect(dsn, autocommit=True) as conn:
        with pytest.raises(psycopg.errors.ReadOnlySqlTransaction):
            conn.execute("CREATE TABLE libephemeris._write_probe (id int)")
        conn.execute("SET default_transaction_read_only=off")
        with pytest.raises(psycopg.errors.InsufficientPrivilege):
            conn.execute(
                "INSERT INTO libephemeris.files VALUES (%s,'probe',1,false)",
                ("0" * 64,),
            )
