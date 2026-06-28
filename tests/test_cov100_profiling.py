# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for libephemeris.profiling.

These tests exercise the profiling utilities deterministically and without
any network access. The heavy ``profile_*`` workloads are driven with the
package-level ``calc_ut``/``houses`` entry points monkeypatched to fast stubs
so the profiling control flow is covered without real ephemeris work.
"""

from __future__ import annotations

import libephemeris as ephem
from libephemeris.profiling import (
    FunctionStats,
    ProfileContext,
    ProfileReport,
    identify_optimization_opportunities,
    profile_house_calculations,
    profile_planetary_calculations,
    run_full_profile,
)


def _make_stats(name: str, calls: int, total_time: float, own_time: float):
    """Build a FunctionStats with explicit fields for report tests."""
    return FunctionStats(
        name=name,
        calls=calls,
        total_time=total_time,
        own_time=own_time,
    )


class TestIdentifyHotPaths:
    """Cover the category-assignment branches in ProfileReport._identify_hot_paths."""

    def test_skyfield_and_core_categories(self):
        """A Skyfield-named and a libephemeris-named function are categorized."""
        report = ProfileReport(
            total_time=10.0,
            function_stats=[
                # name contains "jpllib" -> Skyfield Operations (line 152)
                _make_stats("jpllib.py:_at", calls=10, total_time=1.0, own_time=1.0),
                # name contains "planets" -> Libephemeris Core (line 156)
                _make_stats("planets.py:calc", calls=10, total_time=1.0, own_time=1.0),
            ],
        )

        categories = report._identify_hot_paths()

        assert any(
            "jpllib" in fs.name for fs in categories["Skyfield Operations"]
        )
        assert any(
            "planets" in fs.name for fs in categories["Libephemeris Core"]
        )

    def test_summary_uses_hot_paths(self):
        """summary() walks the hot-path categories and renders percentages."""
        report = ProfileReport(
            operation_name="HotPaths",
            total_time=10.0,
            total_calls=20,
            function_stats=[
                _make_stats("skyfield:_at", calls=10, total_time=5.0, own_time=5.0),
                _make_stats("lunar.py:calc", calls=10, total_time=2.0, own_time=2.0),
            ],
        )

        text = report.summary(top_n=5)

        assert "Skyfield Operations" in text
        assert "Libephemeris Core" in text
        assert "% of total time" in text


class TestPrintStats:
    """Cover the print_stats branches on ProfileReport and ProfileContext."""

    def test_report_print_stats_with_stats(self):
        """ProfileReport.print_stats runs sort+print when stats present (188-190)."""
        with ProfileContext(name="print") as p:
            _ = sum(range(500))

        report = p.get_report()
        assert report.stats is not None

        # pstats writes to the Stats object's own stream (the StringIO created
        # in _generate_report), not to stdout. Exercising the call covers the
        # sort+print branch; the captured buffer holds the rendered table.
        report.print_stats(top_n=3, sort_by="cumulative")
        rendered = report.stats.stream.getvalue()
        assert "function calls" in rendered

    def test_report_print_stats_without_stats(self):
        """ProfileReport.print_stats is a no-op when stats is None."""
        report = ProfileReport(stats=None)
        # Should not raise even though there is nothing to print.
        report.print_stats()

    def test_context_print_stats_with_report(self):
        """ProfileContext.print_stats delegates when a report exists (290-291)."""
        with ProfileContext(name="ctx_print") as p:
            _ = sum(range(500))

        # Delegates to the report's print_stats, which renders into the Stats
        # object's StringIO stream rather than stdout.
        p.print_stats(top_n=3, sort_by="cumulative")
        rendered = p.get_report().stats.stream.getvalue()
        assert "function calls" in rendered

    def test_context_print_stats_without_report(self):
        """ProfileContext.print_stats is a no-op before exit (report is None)."""
        p = ProfileContext(name="no_report")
        # _report is None here; should not raise and should not print.
        p.print_stats()


class TestGenerateReportNaming:
    """Cover the builtin-name else-branch in _generate_report (line 247)."""

    def test_builtin_function_name_uses_raw_funcname(self):
        """Profiling a builtin yields a stat whose name has no file prefix."""
        with ProfileContext(name="builtins") as p:
            for _ in range(5):
                _ = sum(range(200))

        report = p.get_report()
        # Built-in entries have a synthetic filename like "~" / "<built-in>",
        # so their FunctionStats.name is the bare func_name without "file:".
        names = [fs.name for fs in report.function_stats]
        assert any("built-in" in n for n in names)


class TestProfilePlanetaryCalculations:
    """Cover profile_planetary_calculations body (337-382)."""

    def test_default_planets(self, monkeypatch):
        """Default planet list path with speed flag, stubbed calc_ut."""
        calls: list = []

        def fake_calc_ut(jd, planet, flags):
            calls.append((jd, planet, flags))
            return ([0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 0)

        monkeypatch.setattr(ephem, "calc_ut", fake_calc_ut)

        report = profile_planetary_calculations(iterations=2, include_speed=True)

        assert isinstance(report, ProfileReport)
        # 10 default planets x (2 warmup-implicit + 2 iterations) -> calc_ut called.
        assert len(calls) > 0
        # FLG_SPEED should have been passed through.
        from libephemeris.constants import FLG_SPEED

        assert all(flags == FLG_SPEED for _, _, flags in calls)

    def test_explicit_planets_no_speed(self, monkeypatch):
        """Explicit planet subset with include_speed False -> flags == 0."""
        from libephemeris.constants import SUN, MOON

        flags_seen: list = []

        def fake_calc_ut(jd, planet, flags):
            flags_seen.append(flags)
            return ([1.0, 2.0, 3.0], 0)

        monkeypatch.setattr(ephem, "calc_ut", fake_calc_ut)

        report = profile_planetary_calculations(
            iterations=2,
            include_speed=False,
            planets=[SUN, MOON],
        )

        assert isinstance(report, ProfileReport)
        assert flags_seen and all(f == 0 for f in flags_seen)


class TestProfileHouseCalculations:
    """Cover profile_house_calculations body (399-418)."""

    def test_default_house_systems(self, monkeypatch):
        """Default house-system list path, stubbed houses()."""
        seen: list = []

        def fake_houses(jd, lat, lon, hsys):
            seen.append((jd, lat, lon, hsys))
            return ([0.0] * 13, [0.0] * 10)

        monkeypatch.setattr(ephem, "houses", fake_houses)

        report = profile_house_calculations(iterations=2)

        assert isinstance(report, ProfileReport)
        assert len(seen) > 0

    def test_explicit_house_systems(self, monkeypatch):
        """Explicit house-system subset is honored."""
        seen: list = []

        def fake_houses(jd, lat, lon, hsys):
            seen.append(hsys)
            return ([0.0] * 13, [0.0] * 10)

        monkeypatch.setattr(ephem, "houses", fake_houses)

        report = profile_house_calculations(iterations=1, house_systems=["P", "K"])

        assert isinstance(report, ProfileReport)
        # Only P and K (plus the warmup P) should appear.
        assert set(seen) == {ord("P"), ord("K")}


class TestIdentifyOptimizationOpportunities:
    """Cover the suggestion branches in identify_optimization_opportunities."""

    def test_below_threshold_loops_back(self):
        """A function under 10% of time takes the 439->437 fall-through arc."""
        report = ProfileReport(
            total_time=100.0,
            function_stats=[
                # 5% of total -> pct < 10, loop back to next iteration (439->437)
                _make_stats("cool_fn", calls=1, total_time=5.0, own_time=5.0),
            ],
        )

        suggestions = identify_optimization_opportunities(report)

        # No hotspot suggestion for the sub-threshold function; fallback applies.
        assert suggestions == ["No obvious optimization opportunities identified"]

    def test_frame_latlon_hotspot(self):
        """A frame_latlon hotspot hits the dedicated branch (line 446)."""
        # "frame_latlon" contains the substring "at", so the name must NOT also
        # contain "skyfield" or it would match the earlier branch (line 440).
        report = ProfileReport(
            total_time=10.0,
            function_stats=[
                _make_stats(
                    "coordinates:frame_latlon",
                    calls=1,
                    total_time=5.0,
                    own_time=5.0,
                ),
            ],
        )

        suggestions = identify_optimization_opportunities(report)

        assert any("Coordinate transformations" in s for s in suggestions)

    def test_calc_body_hotspot(self):
        """A _calc_body hotspot hits the dedicated branch (line 451)."""
        report = ProfileReport(
            total_time=10.0,
            function_stats=[
                _make_stats(
                    "planets:_calc_body",
                    calls=1,
                    total_time=4.0,
                    own_time=4.0,
                ),
            ],
        )

        suggestions = identify_optimization_opportunities(report)

        assert any("Core calculation function" in s for s in suggestions)

    def test_cache_candidate(self):
        """A frequently-called, slow function is flagged as cache candidate (464)."""
        report = ProfileReport(
            total_time=1000.0,
            function_stats=[
                # >1000 calls and >10us/call, but only 0.2% of total -> not a hotspot.
                _make_stats(
                    "warm_loop",
                    calls=2000,
                    total_time=2.0,
                    own_time=2.0,  # 2.0 / 2000 * 1e6 = 1000 us/call
                ),
            ],
        )

        suggestions = identify_optimization_opportunities(report)

        assert any("CACHE CANDIDATE" in s for s in suggestions)

    def test_io_in_loop(self):
        """An I/O function called many times is flagged (473-474)."""
        report = ProfileReport(
            total_time=1000.0,
            function_stats=[
                # name matches "read" via get_functions_by_name; >10 calls; tiny pct.
                _make_stats(
                    "spk_reader:read_segment",
                    calls=50,
                    total_time=1.0,
                    own_time=1.0,
                ),
            ],
        )

        suggestions = identify_optimization_opportunities(report)

        assert any("I/O IN LOOP" in s for s in suggestions)

    def test_no_opportunities_fallback(self):
        """A benign report yields the no-opportunities fallback (line 480)."""
        report = ProfileReport(
            total_time=100.0,
            function_stats=[
                # <10% pct, <=1000 calls, no I/O name -> no suggestion accumulated.
                _make_stats("tiny", calls=3, total_time=1.0, own_time=1.0),
            ],
        )

        suggestions = identify_optimization_opportunities(report)

        assert suggestions == ["No obvious optimization opportunities identified"]


class TestRunFullProfile:
    """Cover run_full_profile body (496-522)."""

    def _patch_backends(self, monkeypatch):
        """Patch calc_ut and houses to fast deterministic stubs."""

        def fake_calc_ut(jd, planet, flags):
            return ([0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 0)

        def fake_houses(jd, lat, lon, hsys):
            return ([0.0] * 13, [0.0] * 10)

        monkeypatch.setattr(ephem, "calc_ut", fake_calc_ut)
        monkeypatch.setattr(ephem, "houses", fake_houses)

    def test_verbose_prints(self, monkeypatch, capsys):
        """Verbose run prints headers, summaries and suggestions."""
        self._patch_backends(monkeypatch)

        reports = run_full_profile(iterations=2, verbose=True)

        assert set(reports.keys()) == {"planets", "houses"}
        out = capsys.readouterr().out
        assert "Profiling planetary calculations..." in out
        assert "Profiling house calculations..." in out
        assert "COMPREHENSIVE PROFILING RESULTS" in out
        assert "OPTIMIZATION SUGGESTIONS" in out

    def test_quiet(self, monkeypatch, capsys):
        """Non-verbose run produces reports without printing the report block."""
        self._patch_backends(monkeypatch)

        reports = run_full_profile(iterations=1, verbose=False)

        assert set(reports.keys()) == {"planets", "houses"}
        out = capsys.readouterr().out
        assert "COMPREHENSIVE PROFILING RESULTS" not in out
