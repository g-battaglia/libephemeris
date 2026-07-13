"""Tests for libephemeris.tracing module.

Verifies the ContextVar-based accumulator that records which sub-backend
(LEB, Skyfield, Horizons, SPK, ASSIST, Keplerian, Analytical, ERFA, Mixed)
computed each body.
"""

from __future__ import annotations

import logging
import threading

import libephemeris as ephem
import libephemeris.tracing as tracing_module
import pytest
from libephemeris.tracing import (
    _record,
    _record_alias,
    _restore_record,
    _snapshot_record,
    _trace_data,
    get_trace_results,
    start_tracing,
)
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    REGULUS,
    SPICA_STAR,
    FLG_SPEED,
)


class _RecordingHandler(logging.Handler):
    """Collect records attached directly to the installed library logger."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class _CountingContextVar:
    """Proxy that counts ContextVar operations without changing semantics."""

    def __init__(self, inner) -> None:
        self.inner = inner
        self.get_count = 0
        self.set_count = 0
        self.reset_count = 0

    def get(self, *args, **kwargs):
        self.get_count += 1
        return self.inner.get(*args, **kwargs)

    def set(self, *args, **kwargs):
        self.set_count += 1
        return self.inner.set(*args, **kwargs)

    def reset(self, *args, **kwargs):
        self.reset_count += 1
        return self.inner.reset(*args, **kwargs)


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

    def test_record_alias_transfers_internal_source_and_can_restore_previous(self):
        """Aliases expose caller IDs without leaking internal dispatch keys."""
        token = start_tracing()
        try:
            _record(17, "SPK")
            _record_alias(10001, 17)
            _record_alias(10002, 99)
            assert get_trace_results() == {10001: "SPK"}

            _record(17, "LEB")
            previous = _snapshot_record(17)
            _record(17, "ASSIST")
            _record_alias(10002, 17, previous)
            assert get_trace_results() == {
                17: "LEB",
                10001: "SPK",
                10002: "ASSIST",
            }
        finally:
            token.var.reset(token)

        _record_alias(10001, 17)
        assert get_trace_results() == {}

    def test_restore_record_reverts_or_removes_nested_results(self):
        """Failed nested calls can restore both present and absent entries."""
        token = start_tracing()
        try:
            _record(17, "LEB")
            previous = _snapshot_record(17)
            absent = _snapshot_record(18)
            _record(17, "SPK")
            _record(18, "SPK")
            _restore_record(17, previous)
            _restore_record(18, absent)
            assert get_trace_results() == {17: "LEB"}
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

    @pytest.mark.parametrize(
        "entrypoint",
        ["fixstar_ut", "fixstar", "fixstar2_ut", "fixstar2"],
    )
    def test_every_public_fixed_star_entrypoint_records_source(self, entrypoint):
        """All four single-star APIs publish their resolved catalog ID."""
        token = start_tracing()
        try:
            getattr(ephem, entrypoint)("Regulus", self.JD, 0)
            assert get_trace_results()[REGULUS] in ("LEB", "Skyfield", "Mixed")
        finally:
            token.var.reset(token)

    def test_fixed_star_batch_records_each_success(self):
        """The batch API records each successfully returned star."""
        token = start_tracing()
        try:
            ephem.batch_fixstars_ut(
                ("Regulus", "Not A Star", "Spica"), self.JD, 0, skip_errors=True
            )
            traces = get_trace_results()
            assert traces[REGULUS] in ("LEB", "Skyfield")
            assert traces[SPICA_STAR] in ("LEB", "Skyfield")
        finally:
            token.var.reset(token)

    def test_calc_fixed_star_finalization_failure_keeps_previous_trace(
        self, monkeypatch
    ):
        """The outer calc path records only after output finalization succeeds."""
        from libephemeris import planets

        def _fail_finalization(*args, **kwargs):
            raise RuntimeError("forced finalization failure")

        monkeypatch.setattr(planets, "_finalize_output_flags", _fail_finalization)
        token = start_tracing()
        try:
            _record(REGULUS, "STALE")
            with pytest.raises(RuntimeError, match="forced finalization failure"):
                ephem.calc_ut(self.JD, REGULUS, 0)
            assert get_trace_results()[REGULUS] == "STALE"
        finally:
            token.var.reset(token)

    def test_fixed_epoch_transform_failure_keeps_previous_trace(self, monkeypatch):
        """Recursive fixed-frame reduction does not publish partial success."""
        from libephemeris import sidereal_epoch

        def _fail_transform(*args, **kwargs):
            raise RuntimeError("forced transform failure")

        monkeypatch.setattr(
            sidereal_epoch, "transform_fixed_epoch_result", _fail_transform
        )
        ephem.set_sid_mode(ephem.SIDM_J2000)
        token = start_tracing()
        try:
            _record(REGULUS, "STALE")
            with pytest.raises(RuntimeError, match="forced transform failure"):
                ephem.fixstar_ut("Regulus", self.JD, ephem.FLG_SIDEREAL)
            assert get_trace_results()[REGULUS] == "STALE"
        finally:
            token.var.reset(token)
            ephem.reset_session()

    @pytest.mark.parametrize("api_kind", ["module", "context"])
    def test_calc_fixed_epoch_failure_restores_nested_trace(
        self, monkeypatch, api_kind
    ):
        """Module and context recursion restore the caller's prior trace."""
        from libephemeris import sidereal_epoch

        def _fail_transform(*args, **kwargs):
            raise RuntimeError("forced transform failure")

        monkeypatch.setattr(
            sidereal_epoch, "transform_fixed_epoch_result", _fail_transform
        )
        if api_kind == "module":
            ephem.set_sid_mode(ephem.SIDM_J2000)
            api = ephem
        else:
            api = ephem.EphemerisContext()
            api.set_sid_mode(ephem.SIDM_J2000)

        token = start_tracing()
        try:
            _record(SUN, "STALE")
            with pytest.raises(RuntimeError, match="forced transform failure"):
                api.calc_ut(self.JD, SUN, ephem.FLG_SIDEREAL)
            assert get_trace_results()[SUN] == "STALE"
        finally:
            token.var.reset(token)
            ephem.reset_session()

    @pytest.mark.parametrize("api_kind", ["module", "context"])
    @pytest.mark.parametrize("method", ["calc_ut", "calc"])
    def test_fixed_epoch_postprocessing_failure_restores_trace_and_debug(
        self, monkeypatch, api_kind, method
    ):
        """Nested success remains private until every outer conversion succeeds."""
        from libephemeris import sidereal_epoch, state
        from libephemeris.logging_config import get_logger

        monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)
        monkeypatch.setattr(state, "get_horizons_client", lambda: None)
        monkeypatch.setattr(
            sidereal_epoch, "transform_fixed_epoch_result", lambda *args: ()
        )
        logger = get_logger()
        handler = _RecordingHandler()
        previous_level = logger.level
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
        if api_kind == "module":
            ephem.set_sid_mode(ephem.SIDM_J2000)
            api = ephem
        else:
            api = ephem.EphemerisContext()
            api.set_sid_mode(ephem.SIDM_J2000)

        token = start_tracing()
        try:
            _record(SUN, "STALE")
            with pytest.raises(IndexError):
                getattr(api, method)(self.JD, SUN, ephem.FLG_SIDEREAL)
            assert get_trace_results() == {SUN: "STALE"}
            assert not any(
                "source=" in record.getMessage() for record in handler.records
            )
        finally:
            token.var.reset(token)
            logger.removeHandler(handler)
            logger.setLevel(previous_level)
            ephem.reset_session()

    @pytest.mark.parametrize("api_kind", ["module", "context"])
    @pytest.mark.parametrize("method", ["calc_ut", "calc"])
    def test_fixed_epoch_debug_publishes_once_after_outer_success(
        self, monkeypatch, api_kind, method
    ):
        """A successful recursive frame rewrite has one caller-facing source log."""
        from libephemeris import state
        from libephemeris.logging_config import get_logger

        monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)
        monkeypatch.setattr(state, "get_horizons_client", lambda: None)
        logger = get_logger()
        handler = _RecordingHandler()
        previous_level = logger.level
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
        if api_kind == "module":
            ephem.set_sid_mode(ephem.SIDM_J2000)
            api = ephem
        else:
            api = ephem.EphemerisContext()
            api.set_sid_mode(ephem.SIDM_J2000)

        token = start_tracing()
        try:
            getattr(api, method)(self.JD, SUN, ephem.FLG_SIDEREAL)
            assert get_trace_results() == {SUN: "Skyfield"}
            source_messages = [
                record.getMessage()
                for record in handler.records
                if "source=Skyfield" in record.getMessage()
            ]
            assert len(source_messages) == 1
            assert f"body={SUN} " in source_messages[0]
        finally:
            token.var.reset(token)
            logger.removeHandler(handler)
            logger.setLevel(previous_level)
            ephem.reset_session()

    def test_inactive_fallback_contextvar_work_matches_documented_bound(
        self, monkeypatch
    ):
        """Inactive ordinary fallback uses one trace check and one source scope."""
        from libephemeris import planets, state
        from libephemeris.logging_config import get_logger

        trace_var = tracing_module._trace_data
        source_var = planets._dispatch_source
        trace_token = trace_var.set(None)
        source_token = source_var.set(None)
        trace_probe = _CountingContextVar(trace_var)
        source_probe = _CountingContextVar(source_var)
        monkeypatch.setattr(tracing_module, "_trace_data", trace_probe)
        monkeypatch.setattr(planets, "_dispatch_source", source_probe)
        monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)
        monkeypatch.setattr(state, "get_horizons_client", lambda: None)
        monkeypatch.setattr(get_logger(), "level", logging.WARNING)
        monkeypatch.setattr(
            planets,
            "_calc_body",
            lambda t, ipl, iflag: ((1.0, 2.0, 3.0, 0.0, 0.0, 0.0), iflag),
        )

        try:
            ephem.calc_ut(self.JD, SUN, 0)
            assert trace_probe.get_count == 1
            assert (
                source_probe.set_count,
                source_probe.get_count,
                source_probe.reset_count,
            ) == (1, 1, 1)
        finally:
            source_var.reset(source_token)
            trace_var.reset(trace_token)

    def test_inactive_direct_backend_contextvar_work_matches_documented_bound(
        self, monkeypatch
    ):
        """An inactive direct backend checks log capture and trace state once."""
        from libephemeris import fast_calc, planets, state
        from libephemeris.logging_config import get_logger

        sentinel = object()
        trace_var = tracing_module._trace_data
        capture_var = planets._source_log_capture
        trace_token = trace_var.set(None)
        capture_token = capture_var.set(None)
        trace_probe = _CountingContextVar(trace_var)
        capture_probe = _CountingContextVar(capture_var)
        monkeypatch.setattr(tracing_module, "_trace_data", trace_probe)
        monkeypatch.setattr(planets, "_source_log_capture", capture_probe)
        monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
        monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)
        monkeypatch.setattr(state, "get_horizons_client", lambda: None)
        monkeypatch.setattr(get_logger(), "level", logging.WARNING)
        monkeypatch.setattr(
            fast_calc,
            "fast_calc_ut",
            lambda *args, **kwargs: ((1.0, 2.0, 3.0, 0.0, 0.0, 0.0), 0),
        )

        try:
            ephem.calc_ut(self.JD, SUN, 0)
            assert trace_probe.get_count == 1
            assert capture_probe.get_count == 1
            assert capture_probe.set_count == 0
            assert capture_probe.reset_count == 0
        finally:
            capture_var.reset(capture_token)
            trace_var.reset(trace_token)


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
