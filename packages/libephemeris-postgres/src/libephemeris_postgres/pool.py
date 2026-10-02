# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Lazy process-local, read-only PostgreSQL pool."""

from __future__ import annotations

import os
import threading

from libephemeris import CoefficientSourceError
from psycopg.conninfo import conninfo_to_dict
from psycopg_pool import ConnectionPool

from .config import RuntimeConfig

_pool: ConnectionPool | None = None
_config: RuntimeConfig | None = None
_lock = threading.Lock()


def _after_fork() -> None:
    global _pool, _config, _lock
    _pool, _config = None, None
    _lock = threading.Lock()


os.register_at_fork(after_in_child=_after_fork)


def get_pool(config: RuntimeConfig) -> ConnectionPool:
    """Initialize lazily; never touch inherited parent-owned pool threads."""
    global _pool, _config
    with _lock:
        if _pool is not None and _config == config:
            return _pool
        if _pool is not None:
            _pool.close()
            _pool = None
        try:
            # Reject malformed DSNs before pool workers can log them.
            conninfo_to_dict(config.dsn)
            _pool = ConnectionPool(
                config.dsn,
                min_size=0,
                max_size=config.pool_max,
                timeout=config.timeout_seconds,
                kwargs={
                    "autocommit": True,
                    "connect_timeout": max(1, int(config.timeout_seconds)),
                    "options": "-c default_transaction_read_only=on "
                    f"-c statement_timeout={int(config.timeout_seconds * 1000)}",
                },
                open=True,
            )
            _config = config
            return _pool
        except Exception:
            raise CoefficientSourceError(
                "PostgreSQL coefficient source unavailable"
            ) from None


def reset_pool() -> None:
    """Close this process's pool for controlled shutdown or tests."""
    global _pool, _config
    with _lock:
        if _pool is not None:
            _pool.close()
        _pool, _config = None, None
