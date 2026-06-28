# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Lightweight computation tracing for libephemeris.

Allows callers to discover which sub-backend (LEB, Skyfield, Horizons,
SPK, ASSIST, Keplerian) computed each celestial body, with effectively
zero overhead when tracing is not active.

Usage::

    import libephemeris

    token = libephemeris.start_tracing()
    result = libephemeris.calc_ut(jd, body_id, flags)
    traces = libephemeris.get_trace_results()   # {body_id: "LEB", ...}
    token.var.reset(token)
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
