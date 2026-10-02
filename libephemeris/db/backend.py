# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""DB configuration and operation lifetime, isolated from astronomy.

Provenance:
    Project-authored ownership and configuration infrastructure. Only connection
    configuration survives operations; scientific inputs do not. This module
    defines no astronomical model or coefficient.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass
from functools import wraps
from typing import Any, Callable
from uuid import UUID

from .reader import DBReader
from .store import PostgresStore
from ..exceptions import ConfigurationError

_configuration: tuple[str, str] | None = None
_store: PostgresStore | None = None
_store_dsn: str | None = None
_lock = threading.RLock()
# The engine is synchronous and already uses thread-local frame dispatch.
# Compatibility calls borrow one explicit scoped owner through this routing
# slot. It is operation state, never a persistent scientific cache.
_operation_state = threading.local()


def set_db_config(dsn: str | None, dataset_id: str | None = None) -> None:
    """Select an immutable DB dataset or reset to environment/TOML settings.

    Args:
        dsn: PostgreSQL connection string, or None to clear the override.
        dataset_id: Dataset UUID. Required when dsn is supplied.

    Raises:
        ConfigurationError: Configuration is incomplete or UUID is invalid.
    """
    global _configuration
    configuration = None if dsn is None else _validate_config(dsn, dataset_id)
    if dsn is None and dataset_id is not None:
        raise ConfigurationError("Dataset requires a PostgreSQL connection string")
    with _lock:
        close_db()
        _configuration = configuration


def _validate_config(dsn: str, dataset_id: str | None) -> tuple[str, str]:
    """Validate identifiers without displaying connection secrets.

    Args:
        dsn: PostgreSQL connection string.
        dataset_id: Explicit dataset UUID.

    Returns:
        Connection string and canonical UUID string.

    Raises:
        ConfigurationError: Required settings are missing or malformed.
    """
    if not dsn.strip() or not dataset_id:
        raise ConfigurationError(
            "DB mode requires db_url and an explicit db_dataset UUID"
        )
    try:
        canonical_id = str(UUID(dataset_id))
    except (ValueError, AttributeError):
        raise ConfigurationError("DB dataset must be a valid UUID") from None
    return dsn, canonical_id


def get_db_config() -> tuple[str, str]:
    """Resolve setter, environment and TOML configuration in that order.

    Returns:
        Connection string and canonical dataset UUID.

    Raises:
        ConfigurationError: DB mode has not been configured.
    """
    if _configuration is not None:
        return _configuration
    from .._config_toml import get_str

    dsn = os.environ.get("LIBEPHEMERIS_DB_URL") or get_str("db_url") or ""
    dataset = os.environ.get("LIBEPHEMERIS_DB_DATASET") or get_str("db_dataset")
    return _validate_config(dsn, dataset)


def get_db_reader() -> DBReader:
    """Return operation-owned inputs or a detached reader for explicit use.

    Returns:
        A reader whose records must not be placed in process-wide caches.
        Public calculations manage its lifetime through db_operation.
    """
    active = getattr(_operation_state, "reader", None)
    if active is not None:
        return active
    dsn, dataset = get_db_config()
    return DBReader(get_store(dsn), dataset)


def get_store(dsn: str) -> PostgresStore:
    """Share transport resources, not metadata, across dataset versions.

    Args:
        dsn: Configured PostgreSQL connection string.

    Returns:
        Process-owned lazy pool adapter.
    """
    global _store, _store_dsn
    with _lock:
        if _store is None or _store_dsn != dsn:
            close_db()
            _store = PostgresStore(dsn)
            _store_dsn = dsn
        return _store


@dataclass(slots=True)
class DatabaseOperation:
    """Data-only lifetime record for the Python compatibility adapter.

    The caller explicitly invokes begin_operation and finish_operation around
    its owned scope. The numerical kernel neither reads this routing slot nor
    needs a decorator, thread-local, frame-reader object or operation class.

    Attributes:
        reader: Owned DB/routed inputs, or None when borrowing an outer scope.
        failure: First fatal source error, preventing successful partial output.
        previous_reader: Previous Python engine frame source to restore.
        previous_generation: Generation of that frame source.
        previous_has_nutation: Previous frame-source capability.
    """

    reader: Any = None
    failure: Exception | None = None
    previous_reader: object | None = None
    previous_generation: int = -1
    previous_has_nutation: bool = False


def begin_operation(operation: DatabaseOperation) -> None:
    """Acquire inputs only when this is the outermost DB operation.

    Args:
        operation: Explicit lifetime record owned by the caller.

    Raises:
        DBError: Database configuration, connection or metadata fails.
    """
    from ..state import get_calc_mode

    mode = get_calc_mode()
    owner = getattr(_operation_state, "owner", None)
    if owner is not None and owner.failure is not None:
        raise owner.failure
    if mode not in ("db", "routed") or owner is not None:
        return
    from ..fast_calc import _active_local

    operation.previous_reader = getattr(_active_local, "reader", None)
    operation.previous_generation = getattr(_active_local, "gen", -1)
    operation.previous_has_nutation = getattr(_active_local, "has_nutation", False)
    if mode == "routed":
        from ..routing import RoutedReader, get_tier_routes

        operation.reader = RoutedReader(*get_tier_routes())
    else:
        operation.reader = get_db_reader()
    _operation_state.reader = operation.reader
    _operation_state.owner = operation


def finish_operation(operation: DatabaseOperation) -> None:
    """Release owned inputs and restore frame dispatch, even after errors.

    Args:
        operation: Lifetime record passed to begin_operation. Nested borrowers
            do not release their outer owner's inputs. Finishing is idempotent.
    """
    if operation.reader is None:
        return
    from ..fast_calc import _active_local

    _active_local.reader = operation.previous_reader
    _active_local.gen = operation.previous_generation
    _active_local.has_nutation = operation.previous_has_nutation
    _operation_state.reader = None
    _operation_state.owner = None
    operation.reader.close()
    operation.reader = None


def db_operation(function: Callable[..., Any]) -> Callable[..., Any]:
    """Add DB input ownership without changing public calculation signatures.

    Args:
        function: Existing public calculation entry point.

    Returns:
        Signature-preserving wrapper; non-DB backends remain untouched.
    """

    @wraps(function)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        """Execute a calculation within its DB input lifetime.

        Args:
            *args: Positional arguments passed unchanged to the calculation.
            **kwargs: Keyword arguments passed unchanged to the calculation.

        Returns:
            The unmodified calculation result.
        """
        from ..operations import calculation_session

        with calculation_session():
            return function(*args, **kwargs)

    return wrapped


def close_db() -> None:
    """Release connection resources while preserving explicit configuration.

    Call when workers are idle. In-flight reads may fail if their pool is
    closed concurrently; this never permits a fallback to another source.
    """
    global _store, _store_dsn
    with _lock:
        if _store is not None:
            _store.close()
        _store = None
        _store_dsn = None
