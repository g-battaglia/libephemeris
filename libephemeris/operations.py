# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Public synchronous coefficient scopes and fatal-source error propagation.

Provenance:
    Project-authored lifetime infrastructure. No numerical model or persistent
    scientific cache. The compatibility adapter borrows one explicit owner.
"""

from __future__ import annotations

from contextlib import contextmanager
from functools import wraps
import threading
from typing import Any, Callable, Iterator

_resource_lock = threading.RLock()
_active_operations = 0


def _idle_worker_action(function: Callable[..., Any]) -> Callable[..., Any]:
    """Reject resource mutation while any coefficient owner is active.

    Args:
        function: Configuration setter or resource-closing entry point.

    Returns:
        Signature-preserving wrapper holding the lifecycle lock through mutation.
    """

    @wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        from .exceptions import ConfigurationError

        with _resource_lock:
            if _active_operations:
                raise ConfigurationError(
                    "Coefficient resources cannot change during a calculation session"
                )
            return function(*args, **kwargs)

    return guarded


def record_source_failure(error: Exception, *, source: bool = False) -> None:
    """Invalidate the current scope after a fatal storage failure.

    Args:
        error: Original fatal error, retained only until scope exit.
        source: Whether the exception arose at a storage-facing boundary.
    """
    from .db.backend import _operation_state
    from .exceptions import (
        ConfigurationError,
        DBError,
        NetworkSealedError,
        RoutingDataError,
    )
    from .leb_reader import LEBCorruptionError

    if isinstance(
        error, (DBError, NetworkSealedError, LEBCorruptionError, RoutingDataError)
    ) or (source and isinstance(error, ConfigurationError)):
        owner = getattr(_operation_state, "owner", None)
        if owner is not None and owner.failure is None:
            owner.failure = error


def source_guard(function: Callable[..., Any]) -> Callable[..., Any]:
    """Remember storage errors even when an integration catches them.

    Args:
        function: Storage-facing reader method.

    Returns:
        Signature-preserving adapter.
    """

    @wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        from .db.backend import _operation_state

        owner = getattr(_operation_state, "owner", None)
        if owner is not None and owner.failure is not None:
            raise owner.failure
        try:
            return function(*args, **kwargs)
        except Exception as error:
            from .exceptions import LEBCorruptionError, RoutingDataError
            from .state import get_calc_mode

            if isinstance(error, LEBCorruptionError) and get_calc_mode() == "routed":
                fatal = RoutingDataError(
                    "Declared routed LEB data is unavailable or corrupt"
                )
                record_source_failure(fatal)
                raise fatal from error
            record_source_failure(error, source=True)
            raise

    return guarded


@contextmanager
def calculation_session() -> Iterator[None]:
    """Share coefficient inputs across synchronous calls in this thread.

    Nested scopes borrow the outer owner. Inputs and DB metadata are discarded
    on exit; configuration, local file readers and connection pools survive.
    Do not span an async await or reconfigure the backend inside this scope.

    Yields:
        None. Existing calculation signatures remain unchanged.

    Raises:
        DBError: A fatal database error occurred, even if caught inside scope.
        RoutingDataError: A declared file source failed inside the scope.
        ConfigurationError: Explicit backend configuration is unusable.
    """
    from .db.backend import DatabaseOperation, begin_operation, finish_operation

    operation = DatabaseOperation()
    reader = None
    previous_sources = None
    try:
        begin_operation(operation)
        from .db.backend import _operation_state

        reader = getattr(_operation_state, "reader", None)
        if reader is not None and hasattr(reader, "_sources"):
            # Trace each nested call independently, then propagate its actual
            # inputs to its parent. A nested preparation must not erase an
            # outer derivative/search's earlier remote source contribution.
            previous_sources = reader._sources
            reader._sources = set()
        try:
            yield
        except Exception as error:
            record_source_failure(error)
            if operation.failure is not None:
                raise operation.failure from (
                    error if error is not operation.failure else None
                )
            raise
        else:
            if operation.failure is not None:
                raise operation.failure
    finally:
        if reader is not None and previous_sources is not None:
            previous_sources.update(reader._sources)
            reader._sources = previous_sources
        finish_operation(operation)
