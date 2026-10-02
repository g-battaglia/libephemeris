# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Explicit operation record and free-function ownership contract."""

from __future__ import annotations

from libephemeris import fast_calc
from libephemeris.db import backend
from libephemeris.db.reader import DBReader


def test_nested_records_borrow_and_finishing_is_idempotent(db_runtime):
    """Model explicit scoped ownership and a nested borrow without managers.

    Args:
        db_runtime: Native coefficient transport double.
    """
    outer = backend.DatabaseOperation()
    inner = backend.DatabaseOperation()
    previous_reader = getattr(fast_calc._active_local, "reader", None)
    previous_generation = getattr(fast_calc._active_local, "gen", -1)
    backend.begin_operation(outer)
    owned = outer.reader
    assert isinstance(owned, DBReader)
    try:
        backend.begin_operation(inner)
        assert inner.reader is None
        assert backend.get_db_reader() is owned
        owned.eval_body(0, 2451545.0)
        fast_calc._set_active_reader(owned)
        backend.finish_operation(inner)
        assert backend.get_db_reader() is owned
        assert not owned._closed
        assert owned._segments
        assert db_runtime.metadata_calls == 1
    finally:
        backend.finish_operation(outer)
    assert owned._closed
    assert owned._segments == {}
    assert owned._series == {}
    assert getattr(backend._operation_state, "reader", None) is None
    assert fast_calc._active_local.reader is previous_reader
    assert fast_calc._active_local.gen == previous_generation
    backend.finish_operation(outer)
    assert fast_calc._active_local.reader is previous_reader
