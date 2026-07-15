# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Residual branch-coverage tests left after the parallel coverage batches.

Closes the small set of full-suite gaps that the per-module batches missed
(lines that turned out NOT to be covered by other suites).
"""

from __future__ import annotations

import os


from libephemeris import utils
from libephemeris.utils import (
    SPLIT_DEG_KEEP_DEG,
    SPLIT_DEG_KEEP_SIGN,
    SPLIT_DEG_ROUND_MIN,
    split_deg,
)


class TestPlanetaryMoonsResiduals:
    def test_get_moon_coverage_swallows_segment_errors(self, monkeypatch):
        from libephemeris import planetary_moons as pm

        moon_id = pm.PLMOON_OFFSET + 501  # canonical NAIF (Io)

        class _BadSeg:
            start_jd = "not-a-float"  # float() -> ValueError
            end_jd = 0.0

        class _BadKernel:
            class spk:  # noqa: N801 - mimic kernel.spk.segments
                segments = [_BadSeg()]

        monkeypatch.setattr(pm, "_MOON_SPK_BY_BODY", {moon_id: "fake.bsp"})
        monkeypatch.setattr(pm, "_MOON_SPK_KERNELS", {"fake.bsp": _BadKernel()})

        # The except (ValueError/KeyError/AttributeError) arm -> None.
        assert pm.get_moon_coverage(moon_id) is None

    def test_close_moon_kernels_swallows_close_errors(self, monkeypatch):
        from libephemeris import planetary_moons as pm

        class _BadKernel:
            def close(self):
                raise OSError("close failed")

        monkeypatch.setattr(pm, "_MOON_SPK_KERNELS", {"fake.bsp": _BadKernel()})
        monkeypatch.setattr(pm, "_MOON_SPK_BY_BODY", {0: "fake.bsp"})
        pm.close_moon_kernels()  # must not raise; clears the registries
        assert pm._MOON_SPK_KERNELS == {}


class TestUtilsSplitDegResiduals:
    def test_keep_deg_offset_does_not_bump_integer_degree(self):
        # ROUND_MIN offset on a mid-degree value does not change the integer
        # degree, so the KEEP_DEG suppression branch is NOT taken (false arc).
        ideg, imin, isec, secfr, sign = split_deg(
            15.4, roundflag=SPLIT_DEG_ROUND_MIN | SPLIT_DEG_KEEP_DEG
        )
        assert ideg == 15

    def test_keep_sign_offset_does_not_cross_sign(self):
        # Mid-sign value: the rounding offset does not cross the 30 deg sign
        # boundary, so the KEEP_SIGN suppression branch is NOT taken.
        ideg, imin, isec, secfr, sign = split_deg(
            10.4, roundflag=SPLIT_DEG_ROUND_MIN | SPLIT_DEG_KEEP_SIGN
        )
        assert ideg == 10


class TestContextRegisterRelativeFile:
    def test_relative_spk_file_present_in_cwd(self, tmp_path, monkeypatch):
        from libephemeris import state
        from libephemeris.constants import CHIRON
        from libephemeris.context import EphemerisContext

        # Create a relative file that exists in cwd but NOT in the library path,
        # so the "elif not os.path.exists(spk_file)" arc falls through (file
        # exists -> condition False -> proceed).
        monkeypatch.chdir(tmp_path)
        rel = "ctx_rel_probe.bsp"
        (tmp_path / rel).write_bytes(b"\x00")

        empty_lib = tmp_path / "emptylib"
        empty_lib.mkdir()
        monkeypatch.setattr(state, "get_library_path", lambda: str(empty_lib))
        monkeypatch.setattr(state, "_load_spk_kernel", lambda *a, **k: None)
        monkeypatch.setattr(state, "_SPK_KERNELS", {}, raising=False)

        ctx = EphemerisContext()
        # Reaches the relative-path resolution; kernel ends up None -> context
        # registers locally without error.
        ctx.register_spk_body(CHIRON, rel, 2002060)


class TestProfilingResiduals:
    def test_generate_report_builtin_filename(self, monkeypatch):
        from libephemeris import profiling

        class _FakeStats:
            def __init__(self, *a, **k):
                # filename starting with "<" -> the bare-func-name else branch
                self.stats = {("<string>", 1, "evalcode"): (1, 1, 0.1, 0.2, {})}

            def sort_stats(self, *a, **k):
                return self

        monkeypatch.setattr(profiling.pstats, "Stats", _FakeStats)
        with profiling.ProfileContext("probe") as ctx:
            pass
        names = [fs.name for fs in ctx.get_report().function_stats]
        assert "evalcode" in names

    def test_optimization_io_func_below_threshold(self, monkeypatch):
        from libephemeris import profiling

        # An I/O-named function called few times must NOT trigger the I/O-in-loop
        # suggestion (the "fs.calls > 10" false arc).
        # Two I/O funcs: one below the threshold (loop-continue / false arc),
        # one above it (the suggestion-appending true arc). Large total_time
        # keeps pct < 10 so the hotspot loop adds no noise.
        low = profiling.FunctionStats(
            name="reader_fn", calls=3, total_time=1.0, own_time=1.0
        )
        high = profiling.FunctionStats(
            name="reader_hot", calls=50, total_time=1.0, own_time=1.0
        )
        report = profiling.ProfileReport(
            operation_name="probe",
            function_stats=[low, high],
            total_time=1000.0,
            total_calls=53,
        )

        monkeypatch.setattr(
            profiling.ProfileReport,
            "get_functions_by_name",
            lambda self, pat: [low, high] if pat == "read" else [],
        )
        suggestions = profiling.identify_optimization_opportunities(report)
        assert any("I/O IN LOOP: reader_hot" in s for s in suggestions)
        assert not any("I/O IN LOOP: reader_fn" in s for s in suggestions)

    def test_run_full_profile_verbose(self, monkeypatch):
        from libephemeris import profiling

        rep = profiling.ProfileReport(operation_name="x")
        monkeypatch.setattr(
            profiling, "profile_planetary_calculations", lambda *a, **k: rep
        )
        monkeypatch.setattr(
            profiling, "profile_house_calculations", lambda *a, **k: rep
        )
        # verbose=True drives the reports/suggestions print loops.
        out = profiling.run_full_profile(iterations=1, verbose=True)
        assert set(out) == {"planets", "houses"}


def test_azimuth_module_import_smoke():
    # utils import sanity (keeps the module referenced for the pragma'd path).
    assert hasattr(utils, "azalt")
    assert os.path.isabs(os.getcwd())
