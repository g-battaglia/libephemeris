# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Lightweight computation tracing for libephemeris.

Allows callers to discover which sub-backend (LEB, Skyfield, Horizons,
SPK, ASSIST, Keplerian, Analytical, ERFA, Mixed) computed each celestial body,
with minimal overhead when tracing is not active.

Usage::

    import libephemeris

    token = libephemeris.start_tracing()
    result = libephemeris.calc_ut(jd, body_id, flags)
    traces = libephemeris.get_trace_results()   # {body_id: "LEB", ...}
    token.var.reset(token)

Provenance:
    Project-authored ContextVar bookkeeping. It records a label after a
    registered backend succeeds and cannot modify the value or backend choice.
    The module contains no astronomical model or coefficient.
"""

from __future__ import annotations

from contextvars import ContextVar, Token
from typing import Dict, Optional

_trace_data: ContextVar[Optional[Dict[int, str]]] = ContextVar(
    "libephemeris_trace_data", default=None
)


def start_tracing() -> Token[Optional[Dict[int, str]]]:
    """Activate tracing. Returns a token to reset when done."""
    return _trace_data.set({})


def get_trace_results() -> Dict[int, str]:
    """Return ``{body_id: source}`` collected since ``start_tracing()``."""
    return dict(_trace_data.get() or {})


def _record(body_id: int, source: str) -> None:
    """Record a successful computation (internal use only).

    Called at each success dispatch point in ``planets.py`` and
    ``context.py``.  When tracing is inactive the cost is a single
    ``ContextVar.get(None)`` check (~50 ns).
    """
    d = _trace_data.get(None)
    if d is not None:
        d[body_id] = source


def _is_active() -> bool:
    """Return whether programmatic tracing is active in this context."""
    return _trace_data.get(None) is not None


def _snapshot_record(body_id: int) -> tuple[bool, str | None]:
    """Capture a trace entry so an internal alias calculation can restore it."""
    d = _trace_data.get(None)
    if d is None or body_id not in d:
        return False, None
    return True, d[body_id]


def _restore_record(body_id: int, previous: tuple[bool, str | None]) -> None:
    """Restore or remove one entry after a failed nested calculation."""
    d = _trace_data.get(None)
    if d is None:
        return
    existed, source = previous
    if existed and source is not None:
        d[body_id] = source
    else:
        d.pop(body_id, None)


def _record_alias(
    body_id: int,
    source_body_id: int,
    previous_source: tuple[bool, str | None] = (False, None),
) -> None:
    """Transfer an internal body's source to the public caller ID."""
    d = _trace_data.get(None)
    if d is None or source_body_id not in d:
        return
    d[body_id] = d[source_body_id]
    if body_id == source_body_id:
        return
    _restore_record(source_body_id, previous_source)
