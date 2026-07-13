# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for context/module calc entry-point parity."""

from __future__ import annotations

import logging
import re

import pytest

import libephemeris as le
from libephemeris import fast_calc, state
from libephemeris.tracing import _record, _trace_data, get_trace_results, start_tracing


JD = 2451545.0


@pytest.mark.unit
@pytest.mark.parametrize(
    "flags,want",
    [
        (0, 0),
        (le.FLG_SPEED, le.FLG_SPEED),
        (le.FLG_SPEED3, le.FLG_SPEED3),
        (le.FLG_J2000, le.FLG_J2000 | le.FLG_NONUT),
        (
            le.FLG_HELCTR,
            le.FLG_HELCTR | le.FLG_NOGDEFL | le.FLG_NOABERR,
        ),
        (le.FLG_MOSEPH, le.FLG_MOSEPH),
    ],
)
def test_context_calc_tt_retflag_uses_tt_echo_contract(monkeypatch, flags, want):
    """Context calc(), like module calc(), does not inject default SWIEPH."""
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    _, retflag = le.EphemerisContext().calc(JD, le.MARS, flags)

    assert retflag == want


@pytest.mark.unit
@pytest.mark.parametrize(
    "flags,want",
    [
        (0, 0),
        (le.FLG_SPEED3, le.FLG_SPEED3),
        (le.FLG_MOSEPH | le.FLG_SPEED3, le.FLG_MOSEPH | le.FLG_SPEED3),
    ],
)
def test_context_calc_tt_leb_retflag_uses_tt_echo_contract(monkeypatch, flags, want):
    """The TT echo correction also applies to the context LEB fast path."""
    sentinel = object()
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)

    def _fake_fast_calc_tt(reader, tjd, ipl, iflag, **kwargs):
        assert reader is sentinel
        return (0.0, 0.0, 1.0, 0.0, 0.0, 0.0), iflag

    monkeypatch.setattr(fast_calc, "fast_calc_tt", _fake_fast_calc_tt)

    _, retflag = le.EphemerisContext().calc(JD, le.MARS, flags)

    assert retflag == want


@pytest.mark.unit
@pytest.mark.parametrize(
    "method,flags,want",
    [
        ("calc_ut", 0, le.FLG_SWIEPH),
        ("calc", 0, 0),
        ("calc_ut", le.FLG_J2000, le.FLG_SWIEPH | le.FLG_J2000 | le.FLG_NONUT),
        ("calc", le.FLG_J2000, le.FLG_J2000 | le.FLG_NONUT),
        (
            "calc_ut",
            le.FLG_HELCTR,
            le.FLG_SWIEPH | le.FLG_HELCTR | le.FLG_NOGDEFL | le.FLG_NOABERR,
        ),
        (
            "calc_ut",
            le.FLG_TOPOCTR | le.FLG_BARYCTR,
            le.FLG_SWIEPH | le.FLG_TOPOCTR,
        ),
        (
            "calc",
            le.FLG_TOPOCTR | le.FLG_BARYCTR,
            le.FLG_TOPOCTR,
        ),
        (
            "calc_ut",
            le.FLG_SPEED | le.FLG_SPEED3,
            le.FLG_SWIEPH | le.FLG_SPEED,
        ),
    ],
)
def test_context_ecl_nut_retflag_matches_entrypoint_contract(method, flags, want):
    """ECL_NUT keeps center, implied-bit, speed, and TT/UT echo semantics."""
    _, retflag = getattr(le.EphemerisContext(), method)(JD, le.ECL_NUT, flags)

    assert retflag == want


@pytest.mark.unit
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize("body", [le.MARS, le.ECL_NUT])
def test_barycentric_moseph_error_matches_reference_contract(method, body):
    """Module and context entry points reject raw BARYCTR|MOSEPH."""
    flags = le.FLG_BARYCTR | le.FLG_MOSEPH
    message = f"swisseph.{method}: barycentric Moshier positions are not supported."

    for api in (le, le.EphemerisContext()):
        with pytest.raises(le.Error, match=re.escape(message)) as exc_info:
            getattr(api, method)(JD, body, flags)
        assert str(exc_info.value) == message


@pytest.mark.unit
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize(
    "winner",
    [le.FLG_JPLEPH, le.FLG_SWIEPH, le.FLG_TOPOCTR],
)
def test_ecl_nut_barycentric_moseph_honors_priority(method, winner):
    """A higher-priority ephemeris or center suppresses the MOSEPH error."""
    flags = le.FLG_BARYCTR | le.FLG_MOSEPH | winner

    module_result = getattr(le, method)(JD, le.ECL_NUT, flags)
    context_result = getattr(le.EphemerisContext(), method)(JD, le.ECL_NUT, flags)

    assert context_result == module_result


@pytest.mark.unit
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_context_topocentric_earth_skips_leb(monkeypatch, method):
    """LEB must not expose the observer-to-geocentre vector for Earth."""
    sentinel = object()
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    def _unexpected_leb_call(*args, **kwargs):
        pytest.fail("topocentric Earth must bypass the LEB reducer")

    fast_name = "fast_calc_ut" if method == "calc_ut" else "fast_calc_tt"
    monkeypatch.setattr(fast_calc, fast_name, _unexpected_leb_call)

    ctx = le.EphemerisContext()
    ctx.set_topo(12.0, 42.0, 100.0)
    position, retflag = getattr(ctx, method)(JD, le.EARTH, le.FLG_TOPOCTR)

    assert position == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    expected = le.FLG_TOPOCTR | (le.FLG_SWIEPH if method == "calc_ut" else 0)
    assert retflag == expected


@pytest.mark.unit
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize("mode", ["skyfield", "horizons"])
def test_context_forced_backend_mode_skips_context_local_leb(monkeypatch, method, mode):
    """A context-local reader must not defeat a forced non-LEB backend."""
    ctx = le.EphemerisContext()

    def _unexpected_local_reader():
        pytest.fail(f"{mode} mode must not open a context-local LEB reader")

    monkeypatch.setattr(ctx, "get_leb_reader", _unexpected_local_reader)
    monkeypatch.setattr(state, "get_calc_mode", lambda: mode)
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    context_result = getattr(ctx, method)(JD, le.MARS, le.FLG_SPEED)
    module_result = getattr(le, method)(JD, le.MARS, le.FLG_SPEED)

    assert context_result == module_result


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["auto", "leb"])
def test_state_leb_modes_keep_existing_reader(monkeypatch, mode):
    """The horizons gate must not alter auto/leb reader reuse."""
    sentinel = object()
    monkeypatch.setattr(state, "_CALC_MODE", mode)
    monkeypatch.setattr(state, "_LEB_READER", sentinel)

    assert state.get_leb_reader() is sentinel


@pytest.mark.unit
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_horizons_mode_skips_existing_global_leb_for_all_entrypoints(
    monkeypatch, method
):
    """Forced Horizons bypasses an open global reader in module and context APIs."""
    from libephemeris import horizons_backend

    reader = object()
    client = object()
    horizon_position = (123.0, 4.0, 1.5, 0.1, 0.01, 0.001)
    monkeypatch.setattr(state, "_CALC_MODE", "horizons")
    monkeypatch.setattr(state, "_LEB_READER", reader)
    monkeypatch.setattr(state, "get_horizons_client", lambda: client)

    def _unexpected_leb(*args, **kwargs):
        pytest.fail("horizons mode must bypass the existing global LEB reader")

    monkeypatch.setattr(fast_calc, "fast_calc_ut", _unexpected_leb)
    monkeypatch.setattr(fast_calc, "fast_calc_tt", _unexpected_leb)

    def _fake_horizons(actual_client, jd_ut, body, flags, *args):
        assert actual_client is client
        assert body == le.MARS
        return horizon_position, flags

    monkeypatch.setattr(horizons_backend, "horizons_calc_ut", _fake_horizons)

    ctx = le.EphemerisContext()

    def _unexpected_local_reader():
        pytest.fail("horizons mode must not consult a context-local LEB reader")

    monkeypatch.setattr(ctx, "get_leb_reader", _unexpected_local_reader)

    module_result = getattr(le, method)(JD, le.MARS, le.FLG_SPEED)
    context_result = getattr(ctx, method)(JD, le.MARS, le.FLG_SPEED)

    assert module_result[0] == horizon_position
    assert context_result == module_result
    assert state._LEB_READER is reader


@pytest.mark.unit
def test_context_preserves_specific_minor_body_trace(monkeypatch):
    """The context fallback must not overwrite SPK/ASSIST/Keplerian traces."""
    from libephemeris import planets

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    def _fake_minor_body(t, ipl, iflag, ctx):
        _record(ipl, "SPK")
        return (1.0, 2.0, 3.0, 0.0, 0.0, 0.0), iflag

    monkeypatch.setattr(planets, "_calc_body_with_context", _fake_minor_body)

    token = start_tracing()
    try:
        le.EphemerisContext().calc_ut(JD, le.CHIRON, 0)
        assert get_trace_results()[le.CHIRON] == "SPK"
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api", [le, le.EphemerisContext()])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize(
    "body,source",
    [
        (le.MEAN_NODE, "Analytical"),
        (le.MEAN_APOG, "Analytical"),
        (le.CUPIDO, "Analytical"),
        (-le.MEAN_NODE, "Analytical"),
        (le.TRUE_NODE, "Skyfield"),
        (-le.TRUE_NODE, "Skyfield"),
        (le.OSCU_APOG, "Skyfield"),
        (le.INTP_APOG, "Analytical"),
        (le.INTP_PERG, "Analytical"),
    ],
)
def test_analytical_fallback_overwrites_stale_non_minor_trace(
    monkeypatch, api, method, body, source
):
    """Every fallback reports its real source and obeys last-source-wins."""
    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    token = start_tracing()
    try:
        _record(body, "LEB")
        getattr(api, method)(JD, body, 0)
        assert get_trace_results()[body] == source
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api", [le, le.EphemerisContext()])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_ecl_nut_trace_records_its_actual_backend(api, method):
    """The successful ECL_NUT early return participates in tracing."""
    token = start_tracing()
    try:
        _record(le.ECL_NUT, "LEB")
        getattr(api, method)(JD, le.ECL_NUT, 0)
        assert get_trace_results()[le.ECL_NUT] == "ERFA"
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_ast_offset_alias_trace_uses_callers_body_id(monkeypatch, api_kind, method):
    """A built-in asteroid alias exposes the canonical body's source."""
    from libephemeris import planets

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    def _fake_calc_body(t, ipl, iflag, *args):
        assert ipl == le.CERES
        _record(ipl, "SPK")
        return (0.0, 0.0, 1.0, 0.0, 0.0, 0.0), iflag

    if api_kind == "module":
        monkeypatch.setattr(planets, "_calc_body", _fake_calc_body)
        api = le
    else:
        monkeypatch.setattr(planets, "_calc_body_with_context", _fake_calc_body)
        api = le.EphemerisContext()

    requested = le.AST_OFFSET + 1
    token = start_tracing()
    try:
        _record(requested, "STALE")
        getattr(api, method)(JD, requested, 0)
        traces = get_trace_results()
        assert traces[requested] == "SPK"
        assert le.CERES not in traces
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize("backend", ["LEB", "Horizons"])
@pytest.mark.parametrize(
    "requested,dispatched",
    [(le.AST_OFFSET + 1, le.CERES), (-le.MEAN_NODE, le.MEAN_NODE)],
)
def test_alias_trace_on_direct_backends_uses_callers_id(
    monkeypatch, api_kind, method, backend, requested, dispatched
):
    """LEB/Horizons aliases replace stale state under the public request ID."""
    sentinel = object()
    result = ((10.0, 1.0, 1.0, 0.1, 0.0, 0.0), le.FLG_SPEED)

    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(
        state, "get_leb_reader", lambda: sentinel if backend == "LEB" else None
    )
    monkeypatch.setattr(
        state,
        "get_horizons_client",
        lambda: sentinel if backend == "Horizons" else None,
    )

    if backend == "LEB":

        def _fake_fast(reader, tjd, ipl, iflag, **kwargs):
            assert reader is sentinel
            assert ipl == dispatched
            return result

        monkeypatch.setattr(fast_calc, "fast_calc_ut", _fake_fast)
        monkeypatch.setattr(fast_calc, "fast_calc_tt", _fake_fast)
    else:
        from libephemeris import horizons_backend

        def _fake_horizons(client, tjd, ipl, iflag, *args):
            assert client is sentinel
            assert ipl == dispatched
            return result

        monkeypatch.setattr(horizons_backend, "horizons_calc_ut", _fake_horizons)

    api = le if api_kind == "module" else le.EphemerisContext()
    token = start_tracing()
    try:
        _record(requested, "STALE")
        getattr(api, method)(JD, requested, le.FLG_SPEED)
        traces = get_trace_results()
        expected_source = (
            "Analytical"
            if backend == "Horizons" and dispatched == le.MEAN_NODE
            else backend
        )
        assert traces[requested] == expected_source
        assert dispatched not in traces
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
def test_alias_trace_preserves_a_preexisting_canonical_entry(monkeypatch):
    """Internal alias dispatch must not erase an earlier public calculation."""
    from libephemeris import planets

    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    def _fake_calc_body(t, ipl, iflag):
        _record(ipl, "SPK")
        return (0.0, 0.0, 1.0, 0.0, 0.0, 0.0), iflag

    monkeypatch.setattr(planets, "_calc_body", _fake_calc_body)
    requested = le.AST_OFFSET + 1
    token = start_tracing()
    try:
        _record(le.CERES, "LEB")
        le.calc_ut(JD, requested, 0)
        assert get_trace_results() == {le.CERES: "LEB", requested: "SPK"}
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("source", ["LEB", "Skyfield"])
def test_fixed_star_trace_preserves_inner_backend(monkeypatch, source):
    """The outer calc dispatch must not overwrite the fixed-star backend."""
    from libephemeris import fixed_stars, planets

    sentinel = object()
    monkeypatch.setattr(
        state, "get_leb_reader", lambda: sentinel if source == "LEB" else None
    )
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    if source == "LEB":
        monkeypatch.setattr(
            fast_calc,
            "fast_calc_ut",
            lambda *args, **kwargs: (_ for _ in ()).throw(KeyError("star fallback")),
        )
        monkeypatch.setattr(
            fixed_stars,
            "_calc_star_position_leb",
            lambda *args, **kwargs: (120.0, 1.0, 1.0e9),
        )
    else:
        monkeypatch.setattr(
            fixed_stars,
            "_calc_star_position_skyfield",
            lambda *args, **kwargs: (120.0, 1.0, 1.0e9),
        )

    token = start_tracing()
    try:
        _record(le.REGULUS, "STALE")
        planets.calc_ut(JD, le.REGULUS, 0)
        assert get_trace_results()[le.REGULUS] == source
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
def test_failed_outer_finalization_restores_alias_and_internal_trace(
    monkeypatch, api_kind
):
    """A failed public call cannot leak a completed inner backend trace."""
    from libephemeris import planets

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    def _fake_calc_body(t, ipl, iflag, *args):
        assert ipl == le.CERES
        _record(ipl, "SPK")
        return (0.0, 0.0, 1.0, 0.0, 0.0, 0.0), iflag

    def _fail_finalization(*args, **kwargs):
        raise RuntimeError("forced finalization failure")

    monkeypatch.setattr(planets, "_finalize_output_flags", _fail_finalization)
    if api_kind == "module":
        monkeypatch.setattr(planets, "_calc_body", _fake_calc_body)
        api = le
    else:
        monkeypatch.setattr(planets, "_calc_body_with_context", _fake_calc_body)
        api = le.EphemerisContext()

    requested = le.AST_OFFSET + 1
    token = start_tracing()
    try:
        _record(requested, "STALE")
        _record(le.CERES, "LEB")
        with pytest.raises(RuntimeError, match="forced finalization failure"):
            api.calc_ut(JD, requested, 0)
        assert get_trace_results() == {requested: "STALE", le.CERES: "LEB"}
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
def test_failed_type21_path_does_not_record_success(monkeypatch):
    """Observer validation failure must leave the previous trace untouched."""
    from libephemeris import planets, spk

    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)
    monkeypatch.setattr(planets, "get_topo", lambda: None)
    monkeypatch.setattr(spk, "get_spk_type21_target", lambda *args: object())

    token = start_tracing()
    try:
        _record(le.CERES, "STALE")
        with pytest.raises(le.ConfigurationError):
            le.calc_ut(JD, le.CERES, le.FLG_TOPOCTR)
        assert get_trace_results()[le.CERES] == "STALE"
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api", [le, le.EphemerisContext()])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize("body", [le.TRUE_NODE, le.OSCU_APOG])
@pytest.mark.parametrize("flags", [le.FLG_HELCTR, le.FLG_BARYCTR])
def test_centerless_lunar_point_zero_uses_analytical_trace(
    monkeypatch, api, method, body, flags
):
    """A local zero convention must not be attributed to the JPL pipeline."""
    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)

    token = start_tracing()
    try:
        result, _retflag = getattr(api, method)(JD, body, flags)
        assert result == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        assert get_trace_results() == {body: "Analytical"}
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize(
    "body,flags",
    [
        (le.SUN, le.FLG_HELCTR),
        (le.EARTH, 0),
        (le.EARTH, le.FLG_TOPOCTR),
    ],
)
def test_degenerate_origins_report_analytical_source(
    monkeypatch, caplog, api_kind, method, body, flags
):
    """Local zero conventions must not claim a Skyfield vector evaluation."""
    from libephemeris import planets
    from libephemeris.logging_config import get_logger

    class _NoVectorEvaluation:
        def at(self, *args, **kwargs):
            pytest.fail("the degenerate origin must not evaluate a Skyfield vector")

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)
    monkeypatch.setattr(
        planets, "get_planets", lambda: {"earth": _NoVectorEvaluation()}
    )
    monkeypatch.setattr(get_logger(), "propagate", True)

    api = le if api_kind == "module" else le.EphemerisContext()
    if flags & le.FLG_TOPOCTR:
        api.set_topo(0.0, 0.0, 0.0)

    token = start_tracing()
    try:
        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            result, _retflag = getattr(api, method)(JD, body, flags)
        assert result == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        assert get_trace_results() == {body: "Analytical"}
        source_messages = [
            message for message in caplog.messages if "source=Analytical" in message
        ]
        assert len(source_messages) == 1
        assert f"body={body} " in source_messages[0]
        assert not any("source=Skyfield" in message for message in caplog.messages)
    finally:
        _trace_data.reset(token)
        le.reset_session()


@pytest.mark.unit
@pytest.mark.parametrize("api", [le, le.EphemerisContext()])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize(
    "body,flags",
    [
        (le.MEAN_NODE, 0),
        (le.MEAN_APOG, 0),
        (le.CUPIDO, le.FLG_HELCTR),
    ],
)
def test_horizons_local_branches_report_analytical_source(
    monkeypatch, api, method, body, flags
):
    """No-HTTP helpers inside the Horizons dispatcher remain Analytical."""
    sentinel = object()
    monkeypatch.setattr(state, "get_calc_mode", lambda: "horizons")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: sentinel)

    token = start_tracing()
    try:
        getattr(api, method)(JD, body, flags)
        assert get_trace_results() == {body: "Analytical"}
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_failed_direct_leb_output_emits_no_success_source(
    monkeypatch, caplog, api_kind, method
):
    """Neither trace dictionaries nor DEBUG logs publish a failed LEB call."""
    from libephemeris import planets
    from libephemeris.logging_config import get_logger

    sentinel = object()
    result = ((1.0, 2.0, 3.0, 0.0, 0.0, 0.0), 0)
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)
    monkeypatch.setattr(fast_calc, "fast_calc_ut", lambda *args, **kwargs: result)
    monkeypatch.setattr(fast_calc, "fast_calc_tt", lambda *args, **kwargs: result)
    monkeypatch.setattr(
        planets,
        "_to_native_floats",
        lambda values: (_ for _ in ()).throw(RuntimeError("output conversion")),
    )
    monkeypatch.setattr(get_logger(), "propagate", True)
    api = le if api_kind == "module" else le.EphemerisContext()

    token = start_tracing()
    try:
        _record(le.SUN, "STALE")
        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            with pytest.raises(RuntimeError, match="output conversion"):
                getattr(api, method)(JD, le.SUN, 0)
        assert get_trace_results() == {le.SUN: "STALE"}
        assert not any("source=LEB" in message for message in caplog.messages)
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_direct_alias_debug_log_uses_public_body_id(
    monkeypatch, caplog, api_kind, method
):
    """DEBUG output, like structured tracing, exposes the caller's alias ID."""
    from libephemeris.logging_config import get_logger

    sentinel = object()
    result = ((1.0, 2.0, 3.0, 0.0, 0.0, 0.0), 0)
    monkeypatch.setattr(state, "get_calc_mode", lambda: "auto")
    monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)
    monkeypatch.setattr(fast_calc, "fast_calc_ut", lambda *args, **kwargs: result)
    monkeypatch.setattr(fast_calc, "fast_calc_tt", lambda *args, **kwargs: result)
    monkeypatch.setattr(get_logger(), "propagate", True)
    api = le if api_kind == "module" else le.EphemerisContext()
    requested = le.AST_OFFSET + 1

    token = start_tracing()
    try:
        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            getattr(api, method)(JD, requested, 0)
        assert get_trace_results() == {requested: "LEB"}
        source_messages = [m for m in caplog.messages if "source=LEB" in m]
        assert any(f"body={requested} " in message for message in source_messages)
        assert not any(f"body={le.CERES} " in message for message in source_messages)
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
def test_south_node_debug_log_uses_public_id_once(
    monkeypatch, caplog, api_kind, method
):
    """Recursive north-node work is published once under the south-node ID."""
    from libephemeris.logging_config import get_logger

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)
    monkeypatch.setattr(get_logger(), "propagate", True)
    api = le if api_kind == "module" else le.EphemerisContext()
    requested = -le.MEAN_NODE

    token = start_tracing()
    try:
        with caplog.at_level(logging.DEBUG, logger="libephemeris"):
            getattr(api, method)(JD, requested, 0)
        assert get_trace_results() == {requested: "Analytical"}
        source_messages = [m for m in caplog.messages if "source=Analytical" in m]
        assert len(source_messages) == 1
        assert f"body={requested} " in source_messages[0]
    finally:
        _trace_data.reset(token)


@pytest.mark.unit
@pytest.mark.parametrize("api_kind", ["module", "context"])
@pytest.mark.parametrize("method", ["calc_ut", "calc"])
@pytest.mark.parametrize("finalization_fails", [False, True])
def test_keplerian_warning_is_default_level_and_success_only(
    monkeypatch, caplog, api_kind, method, finalization_fails
):
    """Keplerian degradation warns without DEBUG, but never on failed output."""
    from libephemeris import planets
    from libephemeris.logging_config import get_logger

    monkeypatch.setattr(state, "get_calc_mode", lambda: "skyfield")
    monkeypatch.setattr(state, "get_leb_reader", lambda: None)
    monkeypatch.setattr(state, "get_horizons_client", lambda: None)
    monkeypatch.setattr(get_logger(), "propagate", True)

    def _fake_calc_body(t, ipl, iflag, *args):
        planets._mark_dispatch_source("Keplerian")
        return (1.0, 2.0, 3.0, 0.0, 0.0, 0.0), iflag

    if api_kind == "module":
        monkeypatch.setattr(planets, "_calc_body", _fake_calc_body)
        api = le
    else:
        monkeypatch.setattr(planets, "_calc_body_with_context", _fake_calc_body)
        api = le.EphemerisContext()
    if finalization_fails:
        monkeypatch.setattr(
            planets,
            "_finalize_output_flags",
            lambda *args, **kwargs: (_ for _ in ()).throw(
                RuntimeError("output finalization")
            ),
        )

    with caplog.at_level(logging.WARNING, logger="libephemeris"):
        if finalization_fails:
            with pytest.raises(RuntimeError, match="output finalization"):
                getattr(api, method)(JD, le.CHIRON, 0)
        else:
            getattr(api, method)(JD, le.CHIRON, 0)

    warnings = [m for m in caplog.messages if "source=Keplerian" in m]
    assert len(warnings) == (0 if finalization_fails else 1)
