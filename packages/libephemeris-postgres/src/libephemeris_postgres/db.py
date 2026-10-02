# SPDX-License-Identifier: AGPL-3.0-only
from __future__ import annotations

import os
import threading
from typing import Any

import psycopg


class CoefficientSourceError(Exception):
    """Fatal source failure, distinct from scientific range/convergence errors."""


_connection: Any = None
_lock = threading.Lock()


def _after_fork() -> None:
    global _connection, _lock
    if _connection is not None and not _connection.closed:
        # Detach the inherited socket before libpq cleanup can send Terminate.
        os.close(_connection.fileno())
        _connection.close()
    _connection = None
    _lock = threading.Lock()


os.register_at_fork(after_in_child=_after_fork)


def query(sql: str, params: Any = None) -> list[Any]:
    """Serialize reads on one lazy, read-only connection per process."""
    global _connection
    with _lock:
        try:
            if _connection is None or _connection.closed:
                dsn = os.environ.get("LIBEPHEMERIS_PG_URL", "").strip()
                if not dsn:
                    raise ValueError("Missing runtime URL")
                _connection = psycopg.connect(
                    dsn,
                    autocommit=True,
                    connect_timeout=5,
                    options="-c default_transaction_read_only=on -c statement_timeout=5000",
                )
            return _connection.execute(sql, params).fetchall()
        except Exception:
            if _connection is not None:
                _connection.close()
            _connection = None
            raise CoefficientSourceError(
                "PostgreSQL coefficient source unavailable"
            ) from None


def ping() -> None:
    """Probe the transport even when all requested bytes are cached."""
    query("SELECT 1")
