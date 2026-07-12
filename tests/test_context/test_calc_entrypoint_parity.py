# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Regression tests for context/module calc entry-point parity."""

from __future__ import annotations

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
