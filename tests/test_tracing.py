"""Tests for libephemeris.tracing module.

Verifies the ContextVar-based accumulator that records which sub-backend
(LEB, Skyfield, Horizons, SPK, ASSIST, Keplerian) computed each body.
"""

from __future__ import annotations

import threading

import libephemeris as ephem
from libephemeris.tracing import _record, _trace_data, start_tracing, get_trace_results
from libephemeris.constants import SUN, MOON, MERCURY, MARS, FLG_SPEED


# ============================================================================
# Unit tests for the tracing primitives
# ============================================================================


class TestTracingPrimitives:
    """Tests for start_tracing / get_trace_results / _record."""

    def test_no_tracing_by_default(self):
        """When tracing is not started, get_trace_results returns empty dict."""
        # Ensure clean state
        _trace_data.set(None)
        assert get_trace_results() == {}

    def test_start_tracing_returns_token(self):
        """start_tracing returns a contextvars.Token."""
        token = start_tracing()
        try:
            assert token is not None
            assert get_trace_results() == {}
        finally:
            token.var.reset(token)

    def test_record_when_active(self):
        """_record stores body_id -> source when tracing is active."""
        token = start_tracing()
        try:
            _record(0, "LEB")
            _record(1, "Skyfield")
            results = get_trace_results()
            assert results == {0: "LEB", 1: "Skyfield"}
        finally:
            token.var.reset(token)

    def test_record_when_inactive(self):
        """_record is a no-op when tracing is not active."""
        _trace_data.set(None)
        _record(0, "LEB")  # should not raise
        assert get_trace_results() == {}

    def test_record_overwrites(self):
        """Later _record calls overwrite earlier ones for same body_id."""
        token = start_tracing()
        try:
            _record(0, "LEB")
            _record(0, "Skyfield")
            assert get_trace_results()[0] == "Skyfield"
        finally:
            token.var.reset(token)

    def test_token_reset_clears_tracing(self):
        """After token.var.reset(), tracing is inactive again."""
        token = start_tracing()
        _record(0, "LEB")
        token.var.reset(token)
        assert get_trace_results() == {}

    def test_get_trace_results_returns_copy(self):
        """get_trace_results returns a copy, not the internal dict."""
        token = start_tracing()
        try:
            _record(0, "LEB")
            r1 = get_trace_results()
            r1[99] = "fake"
            r2 = get_trace_results()
            assert 99 not in r2
        finally:
            token.var.reset(token)

    def test_nested_tracing(self):
        """Nested start_tracing creates independent sessions."""
        token1 = start_tracing()
        _record(0, "LEB")

        token2 = start_tracing()
        _record(1, "Skyfield")
        inner = get_trace_results()
        assert inner == {1: "Skyfield"}
        token2.var.reset(token2)

        # Outer session is restored
        outer = get_trace_results()
        assert outer == {0: "LEB"}
        token1.var.reset(token1)


# ============================================================================
# Integration tests — tracing through calc_ut
# ============================================================================


class TestTracingIntegration:
    """Verify that calc_ut actually records trace data."""

    # J2000.0
    JD = 2451545.0

    def test_sun_is_traced(self):
        """Computing Sun position should record a trace entry."""
        token = start_tracing()
        try:
            ephem.calc_ut(self.JD, SUN, FLG_SPEED)
            results = get_trace_results()
            assert SUN in results
            assert results[SUN] in ("LEB", "Skyfield", "Horizons")
        finally:
            token.var.reset(token)

    def test_multiple_bodies_traced(self):
        """Computing multiple bodies records all of them."""
        bodies = [SUN, MOON, MERCURY, MARS]
        token = start_tracing()
        try:
            for body in bodies:
                ephem.calc_ut(self.JD, body, FLG_SPEED)
            results = get_trace_results()
            for body in bodies:
                assert body in results, f"body {body} not in trace results"
        finally:
            token.var.reset(token)

    def test_no_trace_when_inactive(self):
        """calc_ut works normally when tracing is off."""
        _trace_data.set(None)
        pos, flags = ephem.calc_ut(self.JD, SUN, FLG_SPEED)
        # Should work fine, no trace data
        assert pos[0] != 0.0
        assert get_trace_results() == {}

    def test_tracing_does_not_affect_result(self):
        """Results are identical with and without tracing."""
        pos_without, _ = ephem.calc_ut(self.JD, SUN, FLG_SPEED)

        token = start_tracing()
        try:
            pos_with, _ = ephem.calc_ut(self.JD, SUN, FLG_SPEED)
        finally:
            token.var.reset(token)

        assert pos_without == pos_with


# ============================================================================
# Thread-safety tests
# ============================================================================


class TestTracingThreadSafety:
    """ContextVar provides per-thread isolation automatically."""

    JD = 2451545.0

    def test_threads_have_independent_traces(self):
        """Each thread gets its own trace accumulator."""
        results_by_thread: dict[str, dict[int, str]] = {}
        errors: list[Exception] = []

        def worker(name: str, body: int):
            try:
                token = start_tracing()
                ephem.calc_ut(self.JD, body, FLG_SPEED)
                results_by_thread[name] = get_trace_results()
                token.var.reset(token)
            except Exception as e:
                errors.append(e)

        t1 = threading.Thread(target=worker, args=("t1", SUN))
        t2 = threading.Thread(target=worker, args=("t2", MOON))
        t1.start()
        t2.start()
        t1.join()
        t2.join()

        assert not errors, f"Thread errors: {errors}"
        # Each thread should have only its own body
        assert SUN in results_by_thread["t1"]
        assert MOON not in results_by_thread["t1"]
        assert MOON in results_by_thread["t2"]
        assert SUN not in results_by_thread["t2"]


# ============================================================================
# Public API tests (imported from libephemeris top-level)
# ============================================================================


class TestPublicAPI:
    """Verify tracing is accessible from the top-level package."""

    def test_start_tracing_importable(self):
        assert hasattr(ephem, "start_tracing")
        assert callable(ephem.start_tracing)

    def test_get_trace_results_importable(self):
        assert hasattr(ephem, "get_trace_results")
        assert callable(ephem.get_trace_results)

    def test_roundtrip_via_public_api(self):
        token = ephem.start_tracing()
        try:
            ephem.calc_ut(2451545.0, SUN, FLG_SPEED)
            results = ephem.get_trace_results()
            assert SUN in results
        finally:
            token.var.reset(token)
