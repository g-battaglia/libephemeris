# SPDX-License-Identifier: AGPL-3.0-only
from __future__ import annotations

from contextlib import contextmanager
from types import SimpleNamespace

import pytest
from libephemeris import CoefficientSourceError

from libephemeris_postgres.config import RuntimeConfig, runtime_config
from libephemeris_postgres.source import BLOCK_SIZE, PgByteSource


@pytest.fixture
def storage(monkeypatch):
    data = bytes(range(256)) * (BLOCK_SIZE // 256 * 3) + b"last"
    blocks = {n: data[n * BLOCK_SIZE : (n + 1) * BLOCK_SIZE] for n in range(4)}
    calls = []

    class Connection:
        def execute(self, sql, params):
            calls.append(params[1])
            return SimpleNamespace(
                fetchall=lambda: [(n, blocks[n]) for n in params[1] if n in blocks]
            )

    @contextmanager
    def connection():
        yield Connection()

    monkeypatch.setattr(
        "libephemeris_postgres.source.get_pool",
        lambda config: SimpleNamespace(connection=connection),
    )
    return data, blocks, calls


def test_read_crosses_blocks_and_keeps_local_payload_with_tiny_cache(storage):
    data, _, calls = storage
    source = PgByteSource(
        "medium_core.leb2", "a" * 64, len(data), RuntimeConfig("unused", cache_blocks=1)
    )
    assert source.read(10, len(data) - 10) == data[10:]
    assert calls == [[0, 1, 2, 3]]
    assert len(source._blocks) == 1
    assert source.read(len(data) - 4, 4) == b"last"
    assert len(calls) == 1
    source.close()
    with pytest.raises(CoefficientSourceError):
        source.read(0, 1)


@pytest.mark.parametrize("offset,size", [(-1, 1), (0, -1), (BLOCK_SIZE * 4, 1)])
def test_invalid_range(storage, offset, size):
    data, _, _ = storage
    source = PgByteSource(
        "medium_core.leb2", "a" * 64, len(data), RuntimeConfig("unused")
    )
    with pytest.raises(CoefficientSourceError):
        source.read(offset, size)


@pytest.mark.parametrize("bad", [None, b"short"])
def test_missing_truncated_block_is_fatal(storage, bad):
    data, blocks, _ = storage
    if bad is None:
        del blocks[1]
    else:
        blocks[1] = bad
    source = PgByteSource(
        "medium_core.leb2", "a" * 64, len(data), RuntimeConfig("unused")
    )
    with pytest.raises(CoefficientSourceError):
        source.read(BLOCK_SIZE, 10)


def test_transport_errors_are_sanitized(storage, monkeypatch):
    data, _, _ = storage

    def failed_pool(config):
        raise RuntimeError("postgresql://user:secret@host/db")

    monkeypatch.setattr("libephemeris_postgres.source.get_pool", failed_pool)
    source = PgByteSource(
        "medium_core.leb2", "a" * 64, len(data), RuntimeConfig("unused")
    )
    with pytest.raises(CoefficientSourceError) as failure:
        source.read(0, 1)
    assert "secret" not in str(failure.value)
    assert failure.value.__suppress_context__


def test_config_hides_dsn_and_rejects_nonfinite_timeout(monkeypatch):
    monkeypatch.setenv("LIBEPHEMERIS_PG_URL", "postgresql://user:secret@host/db")
    assert "secret" not in repr(runtime_config())
    monkeypatch.setenv("LIBEPHEMERIS_PG_TIMEOUT_SECONDS", "nan")
    with pytest.raises(Exception, match="finite"):
        runtime_config()
