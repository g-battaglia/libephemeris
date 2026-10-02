# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Lazy, process-safe PostgreSQL connection pool."""

from __future__ import annotations

import hashlib
import os
import threading
from typing import Any

from .config import RuntimeConfig

_pool: Any = None
_pool_pid: int | None = None
_pool_fingerprint: tuple[str, int, float, int] | None = None
_pool_lock = threading.Lock()


def _after_fork() -> None:
    """Forget inherited pool state without touching parent-owned threads."""

    global _pool, _pool_pid, _pool_fingerprint, _pool_lock
    _pool = None
    _pool_pid = None
    _pool_fingerprint = None
    _pool_lock = threading.Lock()


os.register_at_fork(after_in_child=_after_fork)


def _fingerprint(config: RuntimeConfig) -> tuple[str, int, float, int]:
    # Keep only a digest of the DSN: no credentials are retained in diagnostics.
    digest = hashlib.sha256(config.dsn.encode("utf-8")).hexdigest()
    return (digest, config.pool_max, config.timeout_seconds, config.cache_pages)


def _safe_unavailable() -> Exception:
    from libephemeris import CoefficientSourceError

    return CoefficientSourceError("PostgreSQL coefficient source unavailable")


def get_pool(config: RuntimeConfig) -> Any:
    """Create a pool lazily and never reuse one across a fork."""

    global _pool, _pool_pid, _pool_fingerprint
    pid = os.getpid()
    fingerprint = _fingerprint(config)
    with _pool_lock:
        if _pool is not None and _pool_pid == pid and _pool_fingerprint == fingerprint:
            return _pool
        if _pool is not None and _pool_pid == pid and _pool_fingerprint != fingerprint:
            try:
                _pool.close()
            except Exception:
                pass
            _pool = None
            _pool_fingerprint = None
        try:
            from psycopg.conninfo import conninfo_to_dict
            from psycopg_pool import ConnectionPool

            # Parse before constructing the pool so malformed connection info is
            # not logged by psycopg_pool with credentials attached.
            conninfo_to_dict(config.dsn)
            timeout_ms = int(config.timeout_seconds * 1000)
            _pool = ConnectionPool(
                conninfo=config.dsn,
                min_size=0,
                max_size=config.pool_max,
                timeout=config.timeout_seconds,
                kwargs={
                    "autocommit": True,
                    "connect_timeout": max(1, int(config.timeout_seconds)),
                    "options": (
                        "-c default_transaction_read_only=on "
                        f"-c statement_timeout={timeout_ms}"
                    ),
                },
                open=True,
            )
            _pool_pid = pid
            _pool_fingerprint = fingerprint
            return _pool
        except Exception:
            raise _safe_unavailable() from None


def reset_pool() -> None:
    """Close the current process pool; intended for tests and shutdown."""

    global _pool, _pool_pid, _pool_fingerprint
    with _pool_lock:
        if _pool is not None:
            try:
                _pool.close()
            except Exception:
                pass
        _pool = None
        _pool_pid = None
        _pool_fingerprint = None
