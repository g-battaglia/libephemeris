"""Native trace lifecycle, snapshots and thread ownership."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
import inspect
from pathlib import Path
import threading

import libephemeris as eph
from libephemeris import tracing
import pytest


@pytest.fixture(autouse=True)
def configured_engine():
    eph.set_precision_tier("base")
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    token = tracing.start_tracing()
    assert token.reset()
    yield
    tracing.start_tracing().reset()


def test_public_signatures_and_native_token():
    assert list(inspect.signature(tracing.start_tracing).parameters) == []
    assert list(inspect.signature(tracing.get_trace_results).parameters) == []
    assert tracing.start_tracing is eph.start_tracing
    assert tracing.get_trace_results is eph.get_trace_results
    token = tracing.start_tracing()
    assert isinstance(token, tracing.TraceToken)
    assert token.reset() is True
    assert token.reset() is False


def test_records_real_backend_and_returns_owned_snapshots():
    assert tracing.get_trace_results() == {}
    token = tracing.start_tracing()
    eph.calc_ut(2451545.0, eph.MARS)
    snapshot = tracing.get_trace_results()
    assert snapshot[eph.MARS] == "LEB"
    snapshot[9999] = "caller-owned"
    assert 9999 not in tracing.get_trace_results()
    eph.calc(2451545.0, eph.MOON)
    assert tracing.get_trace_results()[eph.MOON] == "LEB"
    assert token.reset()
    assert tracing.get_trace_results() == {}


def test_nested_scope_restores_parent_and_its_still_live_token():
    first = tracing.start_tracing()
    eph.calc_ut(2451545.0, eph.MARS)
    second = tracing.start_tracing()
    assert tracing.get_trace_results() == {}
    assert first.reset() is False
    eph.calc_ut(2451545.0, eph.VENUS)
    assert second.reset() is True
    assert tracing.get_trace_results() == {eph.MARS: "LEB"}
    assert second.reset() is False
    assert first.reset() is True
    assert tracing.get_trace_results() == {}


def test_multiple_nested_scopes_restore_provenance_even_after_an_exception():
    outer = tracing.start_tracing()
    eph.calc_ut(2451545.0, eph.SUN)
    inner = tracing.start_tracing()
    eph.calc_ut(2451545.0, eph.MOON)
    deepest = tracing.start_tracing()
    try:
        eph.calc_ut(2451545.0, eph.MERCURY)
        assert not outer.reset()
        assert not inner.reset()
        raise RuntimeError("caller failure")
    except RuntimeError:
        assert deepest.reset()
    assert tracing.get_trace_results() == {eph.MOON: "LEB"}
    assert inner.reset()
    assert tracing.get_trace_results() == {eph.SUN: "LEB"}
    assert outer.reset()
    assert tracing.get_trace_results() == {}


def test_foreign_token_cannot_reset_worker_and_threads_are_independent():
    parent_token = tracing.start_tracing()
    eph.calc_ut(2451545.0, eph.MARS)
    barrier = threading.Barrier(8)

    def worker(body):
        token = tracing.start_tracing()
        barrier.wait()
        eph.calc_ut(2451545.0, body)
        snapshot = tracing.get_trace_results()
        assert snapshot[body] == "LEB"
        assert parent_token.reset() is False
        assert token.reset() is True
        assert tracing.get_trace_results() == {}
        return snapshot

    bodies = [eph.SUN, eph.MOON, eph.MERCURY, eph.VENUS, eph.MARS,
              eph.JUPITER, eph.SATURN, eph.URANUS]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(worker, bodies))
    assert all(snapshot[body] == "LEB" for snapshot, body in zip(results, bodies))
    assert tracing.get_trace_results()[eph.MARS] == "LEB"
    assert parent_token.reset() is True
