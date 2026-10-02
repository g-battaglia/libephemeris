# SPDX-License-Identifier: AGPL-3.0-only
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import pytest

from libephemeris_postgres import db
from libephemeris_postgres.source import (
    BLOCK_SIZE,
    CoefficientSourceError,
    PgByteSource,
)


@pytest.mark.parametrize("configured", [" MEDIUM, extended, ", "medium,,extended"])
def test_remote_tier_selection_ignores_empty_items(monkeypatch, configured):
    import libephemeris_postgres as provider
    from pathlib import Path

    core = Path(__file__).resolve().parents[3] / "libephemeris/data/leb2/base_core.leb2"
    monkeypatch.setenv("LIBEPHEMERIS_PG_TIERS", configured)
    monkeypatch.setattr(provider, "get_precision_tier", lambda: "base")
    monkeypatch.setattr(
        provider, "_discover_reviewed_leb_tier_cores", lambda: {"base": str(core)}
    )
    with provider.open_reader() as reader:
        assert reader.eval_body(0, 2451545.0)


@pytest.fixture
def storage(monkeypatch):
    data = bytes(range(256)) * (BLOCK_SIZE // 256 * 3) + b"last"
    blocks = {n: data[n * BLOCK_SIZE : (n + 1) * BLOCK_SIZE] for n in range(4)}
    calls = []

    def query(sql, params):
        calls.append(params[1])
        return [(n, blocks[n]) for n in params[1] if n in blocks]

    monkeypatch.setattr("libephemeris_postgres.source.query", query)
    return data, blocks, calls


def test_slices_cache_concurrent_reads_and_close(storage, monkeypatch):
    data, _, calls = storage
    monkeypatch.setattr("libephemeris_postgres.source.CACHE_BLOCKS", 1)
    source = PgByteSource("a" * 64, len(data))
    assert source[10 : len(data)] == data[10:]
    assert calls == [[0, 1, 2, 3]]
    assert len(source._blocks) == 1
    with ThreadPoolExecutor(8) as workers:
        assert (
            list(workers.map(lambda _: source[len(data) - 4 : len(data)], range(32)))
            == [b"last"] * 32
        )
    assert len(calls) == 1
    source.close()
    with pytest.raises(CoefficientSourceError):
        source[0:1]


@pytest.mark.parametrize("start,stop", [(-1, 1), (1, 0), (0, BLOCK_SIZE * 4)])
def test_invalid_range(storage, start, stop):
    data, _, _ = storage
    with pytest.raises(CoefficientSourceError):
        PgByteSource("a" * 64, len(data))[start:stop]


@pytest.mark.parametrize("bad", [None, b"short"])
def test_missing_truncated_block_is_fatal(storage, bad):
    data, blocks, _ = storage
    if bad is None:
        del blocks[1]
    else:
        blocks[1] = bad
    with pytest.raises(CoefficientSourceError):
        PgByteSource("a" * 64, len(data))[BLOCK_SIZE : BLOCK_SIZE + 10]


def test_remote_chunk_corruption_and_range_errors_remain_distinct(monkeypatch):
    import math
    from pathlib import Path
    from libephemeris.leb2_reader import LEB2Reader
    from libephemeris_postgres.source import PgLEB2Reader

    payload = bytearray(
        (
            Path(__file__).resolve().parents[3]
            / "libephemeris/data/leb2/base_core.leb2"
        ).read_bytes()
    )
    with LEB2Reader(
        str(
            Path(__file__).resolve().parents[3]
            / "libephemeris/data/leb2/base_core.leb2"
        )
    ) as local:
        chunk = local._chunk_index[1][0]
        payload[chunk.blob_offset : chunk.blob_offset + chunk.compressed_size] = (
            b"x" * chunk.compressed_size
        )
    monkeypatch.setattr(
        "libephemeris_postgres.source.query",
        lambda sql, params: [
            (key, bytes(payload[key * BLOCK_SIZE : (key + 1) * BLOCK_SIZE]))
            for key in params[1]
        ],
    )
    with PgLEB2Reader(
        "postgres://pin/base_core.leb2", data=PgByteSource("pin", len(payload))
    ) as reader:
        with pytest.raises(ValueError):
            reader.eval_body(1, math.nextafter(reader._bodies[1].jd_start, -math.inf))
        with pytest.raises(CoefficientSourceError):
            reader.eval_body(1, reader._bodies[1].jd_start)
        reader._mm._blocks.clear()

        def unavailable(sql, params):
            raise CoefficientSourceError("PostgreSQL coefficient source unavailable")

        monkeypatch.setattr("libephemeris_postgres.source.query", unavailable)
        with pytest.raises(CoefficientSourceError):
            reader.eval_nutation(reader._nutation.jd_start)


def test_transport_sanitized_reconnect_and_fork(monkeypatch):
    monkeypatch.setattr(db, "_connection", None)
    monkeypatch.setenv("LIBEPHEMERIS_PG_URL", "postgresql://user:secret@host/db")
    calls = []

    class Connection:
        closed = False
        failed = False

        def execute(self, sql, params):
            if self.failed:
                raise RuntimeError("secret")
            return self

        def fetchall(self):
            return [(1,)]

        def close(self):
            self.closed = True

    def connect(*args, **kwargs):
        calls.append(Connection())
        return calls[-1]

    monkeypatch.setattr(db.psycopg, "connect", connect)
    db.ping()
    calls[-1].failed = True
    with pytest.raises(CoefficientSourceError) as failure:
        db.ping()
    assert "secret" not in str(failure.value)
    assert calls[0].closed
    db.ping()
    assert len(calls) == 2
    monkeypatch.setattr(db.os, "close", lambda fd: None)
    monkeypatch.setattr(Connection, "fileno", lambda self: 123, raising=False)
    db._after_fork()
    db.ping()
    assert len(calls) == 3
    assert calls[1].closed
    db._connection.close()
    db._connection = None
