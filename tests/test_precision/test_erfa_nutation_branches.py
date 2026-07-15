# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for libephemeris.erfa_nutation.

pyerfa is installed in the test environment, so the ``_HAS_ERFA is False``
fallback branches never run normally. These tests flip the flag (and, once,
reload the module without erfa) to exercise every path.
"""

from __future__ import annotations

import importlib.util
import sys

import pytest

import libephemeris.erfa_nutation as en

_JD = 2451545.0


class TestErfaAbsentBranches:
    """Runtime fallback branches taken when pyerfa is unavailable."""

    @pytest.fixture(autouse=True)
    def _force_absent(self, monkeypatch):
        monkeypatch.setattr(en, "_HAS_ERFA", False)
        en.clear_erfa_nutation_cache()
        yield
        en.clear_erfa_nutation_cache()

    def test_has_erfa_reports_false(self):
        assert en.has_erfa() is False

    def test_nut00a_returns_none(self):
        assert en.get_erfa_nutation_nut00a(_JD) is None

    def test_nut06a_returns_none(self):
        assert en.get_erfa_nutation_nut06a(_JD) is None

    def test_pnm06a_returns_none(self):
        assert en.get_erfa_pnm06a_matrix(_JD) is None

    def test_obliquity_returns_none(self):
        assert en.get_erfa_obliquity_iau2006(_JD) is None

    def test_cached_falls_back_to_skyfield(self):
        dpsi, deps = en.get_erfa_nutation_cached(_JD + 0.5)
        assert isinstance(dpsi, float)
        assert isinstance(deps, float)
        assert abs(dpsi) < 1e-3 and abs(deps) < 1e-3

    def test_compare_models_skips_erfa_block(self):
        result = en.compare_nutation_models(_JD)
        assert "skyfield_iau2000a" in result
        assert "erfa_nut00a" not in result
        assert "differences_mas" not in result


class TestErfaPresentDefensive:
    """Defensive branches reachable only when an erfa getter returns falsy."""

    def test_cached_raises_when_model_yields_none(self, monkeypatch):
        assert en._HAS_ERFA is True
        en.clear_erfa_nutation_cache()
        monkeypatch.setattr(en, "get_erfa_nutation_nut00a", lambda _jd: None)
        with pytest.raises(ImportError):
            en.get_erfa_nutation_cached(_JD + 0.321, "nut00a")
        en.clear_erfa_nutation_cache()

    def test_compare_models_handles_falsy_getter_results(self, monkeypatch):
        assert en._HAS_ERFA is True
        monkeypatch.setattr(en, "get_erfa_nutation_nut00a", lambda _jd: None)
        monkeypatch.setattr(en, "get_erfa_nutation_nut06a", lambda _jd: None)
        result = en.compare_nutation_models(_JD)
        assert "erfa_nut00a" not in result
        assert "erfa_nut06a" not in result
        assert "differences_mas" not in result


def test_import_time_fallback_when_erfa_missing(monkeypatch):
    """Exec the module source with erfa absent to hit the import except arm.

    Executing into a throwaway module object (rather than reloading the live
    ``libephemeris.erfa_nutation``) still records line coverage for the file
    while leaving the imported module — and any references other tests hold —
    completely untouched.
    """
    monkeypatch.setitem(sys.modules, "erfa", None)
    spec = importlib.util.spec_from_file_location(
        "erfa_nutation_noerfa_probe", en.__file__
    )
    fresh = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fresh)

    assert fresh._HAS_ERFA is False
    assert fresh.has_erfa() is False
    # The live module is unaffected.
    assert en._HAS_ERFA is True
