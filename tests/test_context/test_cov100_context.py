# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.context`.

These tests target the specific uncovered lines and branch arcs in
``EphemerisContext`` that are not exercised by the comprehensive suite:
the LEB-reader close error path, the ``get_sid_mode`` overload stubs and
``full=True`` branch, the SPK registration edge cases, and the
Terrestrial-Time ``calc()`` LEB fast-path fallback plus ``calc_pctr()``.

Inputs are deterministic (fixed Julian Days, Rome coordinates) and all
LEB/HTTP seams are monkeypatched so the file passes with no network.
"""

from __future__ import annotations

import typing

import pytest

from libephemeris import state
from libephemeris.constants import SUN, MOON, MARS
from libephemeris.context import EphemerisContext


# Deterministic test inputs.
J2000 = 2451545.0
ROME_LON = 12.5
ROME_LAT = 41.9


class TestSetLebFileClose:
    """Cover the ``set_leb_file`` reader-close error handling (207-208)."""

    @pytest.mark.unit
    def test_close_swallows_attribute_and_os_error(self):
        """A reader whose ``close`` raises is silently ignored on reset."""

        class _BadReader:
            def close(self):
                raise OSError("boom")

        ctx = EphemerisContext()
        ctx._leb_reader = _BadReader()
        ctx._leb_file = "/tmp/does-not-matter.leb"

        # Must not raise even though close() throws OSError (lines 207-208).
        ctx.set_leb_file(None)

        assert ctx._leb_file is None
        assert ctx._leb_reader is None


class TestGetSidModeOverloads:
    """Cover the ``get_sid_mode`` overload stubs and ``full`` branch."""

    @pytest.mark.unit
    def test_overload_stub_bodies_execute(self):
        """Calling the registered overload stubs runs their ``...`` bodies.

        This exercises the ``289->exit`` and ``292->exit`` branch arcs:
        the ``@overload`` stub function bodies which are never invoked at
        runtime through the public name (the real implementation shadows
        them). We retrieve them via ``typing.get_overloads`` and call
        them directly so the ``...`` body falls through to function exit.
        """
        overloads = typing.get_overloads(EphemerisContext.get_sid_mode)
        assert len(overloads) == 2
        for stub in overloads:
            # The body is a bare ``...`` which returns None.
            assert stub(None, True) is None

    @pytest.mark.unit
    def test_get_sid_mode_full_returns_tuple(self):
        """``get_sid_mode(full=True)`` returns the full config tuple (305)."""
        ctx = EphemerisContext()
        ctx.set_sid_mode(2, t0=2440000.0, ayan_t0=24.5)

        result = ctx.get_sid_mode(full=True)

        assert result == (2, 2440000.0, 24.5)


class TestRegisterSpkBody:
    """Cover SPK registration edge cases (380-381, 391->407, 399->401)."""

    @pytest.mark.unit
    def test_relative_file_not_found_raises(self, monkeypatch):
        """A relative path missing from lib dir and cwd raises (380-381).

        ``full_path`` (lib_path + spk_file) does not exist and the bare
        relative ``spk_file`` does not exist either, so the ``elif`` at
        line 380 is True and line 381 raises.
        """
        monkeypatch.setattr(state, "get_library_path", lambda: "/nonexistent-dir")

        ctx = EphemerisContext()
        with pytest.raises(FileNotFoundError, match="SPK file not found"):
            ctx.register_spk_body(99001, "no_such_file.bsp", 2099001)

    @pytest.mark.unit
    def test_kernel_none_skips_validation(self, monkeypatch, tmp_path):
        """When the cached kernel is None, validation is skipped (391->407).

        The file exists and ``_load_spk_kernel`` is stubbed so nothing is
        stored in ``_SPK_KERNELS``; ``kernel`` is then None and execution
        jumps straight to the context-local registration at line 407.
        """
        spk_path = tmp_path / "fake.bsp"
        spk_path.write_bytes(b"not a real kernel")

        monkeypatch.setattr(state, "_load_spk_kernel", lambda fp: None)
        # Ensure no stale kernel is cached for this path.
        monkeypatch.setitem(state._SPK_KERNELS, str(spk_path), None)
        state._SPK_KERNELS.pop(str(spk_path), None)

        ctx = EphemerisContext()
        ctx.register_spk_body(99002, str(spk_path), 2099002)

        assert ctx.get_spk_body_info(99002) == (str(spk_path), 2099002)

    @pytest.mark.unit
    def test_kernel_without_names_raises_value_error(self, monkeypatch, tmp_path):
        """A kernel missing the body and lacking ``names`` raises (399->401).

        Both ``kernel[naif_id]`` and ``kernel[str(naif_id)]`` raise
        KeyError, and the kernel has no ``names`` attribute, so the
        ``hasattr`` check at line 399 is False and line 401 raises.
        """
        spk_path = tmp_path / "fake2.bsp"
        spk_path.write_bytes(b"not a real kernel")

        class _NamelessKernel:
            def __getitem__(self, key):
                raise KeyError(key)

        kernel = _NamelessKernel()
        assert not hasattr(kernel, "names")

        monkeypatch.setattr(state, "_load_spk_kernel", lambda fp: kernel)
        monkeypatch.setitem(state._SPK_KERNELS, str(spk_path), kernel)

        ctx = EphemerisContext()
        with pytest.raises(ValueError, match="not found in SPK kernel"):
            ctx.register_spk_body(99003, str(spk_path), 2099003)


class TestCalcTtLebPath:
    """Cover the Terrestrial-Time ``calc()`` LEB fast path (565-612)."""

    @pytest.mark.unit
    def test_calc_falls_back_to_global_reader(self, monkeypatch):
        """No context reader -> consult the global reader (567-569).

        With both the context reader and the global reader returning None,
        ``calc()`` skips the LEB block and uses the Skyfield fallback.
        """
        ctx = EphemerisContext()
        ctx._leb_file = None
        ctx._leb_reader = None

        # Force the global reader to None deterministically (the line
        # 567-569 branch that consults state.get_leb_reader()).
        monkeypatch.setattr(state, "get_leb_reader", lambda: None)

        result, retflag = ctx.calc(J2000, SUN, 0)

        assert len(result) == 6
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_calc_leb_topo_then_fallback_on_keyerror(self, monkeypatch):
        """LEB path with topo set, raising KeyError -> fallback (578, 598-612).

        The global reader is a non-None sentinel so the LEB block runs;
        ``self.topo`` is set so the topo tuple is built (line 578); the
        patched ``fast_calc_tt`` raises KeyError, exercising the except
        block (598-605) and the Skyfield fallback (608-612).
        """
        from libephemeris import fast_calc

        ctx = EphemerisContext()
        ctx._leb_file = None
        ctx._leb_reader = None
        ctx.set_topo(ROME_LON, ROME_LAT, 0.0)

        sentinel = object()
        monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)

        def _raise(*args, **kwargs):
            raise KeyError("body not in leb")

        monkeypatch.setattr(fast_calc, "fast_calc_tt", _raise)

        result, retflag = ctx.calc(J2000, MARS, 0)

        assert len(result) == 6
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_calc_leb_topo_then_fallback_on_valueerror(self, monkeypatch):
        """ValueError from the LEB path also triggers the fallback (598-612)."""
        from libephemeris import fast_calc

        ctx = EphemerisContext()
        ctx._leb_file = None
        ctx._leb_reader = None
        ctx.set_topo(ROME_LON, ROME_LAT, 0.0)

        sentinel = object()
        monkeypatch.setattr(state, "get_leb_reader", lambda: sentinel)

        def _raise(*args, **kwargs):
            raise ValueError("jd out of range")

        monkeypatch.setattr(fast_calc, "fast_calc_tt", _raise)

        result, _ = ctx.calc(J2000, MOON, 0)

        assert len(result) == 6
        assert 0 <= result[0] < 360


class TestCalcPctr:
    """Cover the planet-centric calculation method (661-665)."""

    @pytest.mark.unit
    def test_calc_pctr_returns_valid_position(self):
        """``calc_pctr`` delegates to the planets pipeline (661-665)."""
        ctx = EphemerisContext()

        result, retflag = ctx.calc_pctr(J2000, MOON, MARS, 0)

        assert len(result) == 6
        assert 0 <= result[0] < 360
