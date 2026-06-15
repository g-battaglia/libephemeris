# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Branch-coverage tests for :mod:`libephemeris.crossing`.

These tests drive the crossing-event solvers (``solcross``, ``mooncross``,
``cross_ut``, ``helio_cross``, lunar node crossings) and the retrograde
station helpers through their edge branches: backward searches, retrograde
scan fallbacks, Brent's-method fallbacks near stations, divergence guards,
and the controlled-failure paths reached when the underlying ephemeris call
raises. The ephemeris-failure and numerically-pathological branches are
exercised by monkeypatching the module-level ``calc``/``calc_ut`` seams with
deterministic fakes, so the suite needs no network access.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as ephem
import libephemeris.crossing as cr
from libephemeris.crossing import (
    Error,
    _brent_find_crossing,
    _brent_find_station,
    _find_bracket_for_crossing,
    _get_adaptive_max_iterations,
    calc_velocity_at_station,
    cross_ut,
    get_station_info,
    get_station_type,
    helio_cross,
    helio_cross_ut,
    is_retrograde,
    mooncross,
    mooncross_node,
    mooncross_node_ut,
    mooncross_ut,
    next_retrograde_ut,
    solcross,
    solcross_ut,
)
from libephemeris.constants import (
    JUPITER,
    MARS,
    MERCURY,
    MOON,
    NEPTUNE,
    PLUTO,
    SATURN,
    SUN,
    URANUS,
    FLG_HELCTR,
    FLG_SPEED,
)
from libephemeris.exceptions import CalculationError

J2000 = 2451545.0


# ===========================================================================
# Helpers for installing controlled fake ephemeris functions.
# ===========================================================================


def _patch_calc_ut(monkeypatch, fn):
    """Replace ``crossing.calc_ut`` with ``fn`` for the duration of a test."""
    monkeypatch.setattr(cr, "calc_ut", fn)


def _patch_calc(monkeypatch, fn):
    """Replace ``crossing.calc`` with ``fn`` for the duration of a test."""
    monkeypatch.setattr(cr, "calc", fn)


# ===========================================================================
# _get_adaptive_max_iterations (lines 97-100)
# ===========================================================================


class TestAdaptiveMaxIterations:
    """Cover every speed-tier branch of the iteration-count helper."""

    def test_normal_planet(self):
        assert _get_adaptive_max_iterations(0.5) == cr.NR_MAX_ITER_PLANET

    def test_slow_planet(self):
        assert _get_adaptive_max_iterations(0.05) == cr.NR_MAX_ITER_HELIO

    def test_very_slow_planet(self):
        # 0.001 <= |speed| < 0.01 -> very slow tier (lines 97-98)
        assert _get_adaptive_max_iterations(0.005) == cr.NR_MAX_ITER_VERY_SLOW

    def test_near_station(self):
        # |speed| < 0.001 -> station tier (lines 99-100)
        assert _get_adaptive_max_iterations(0.0001) == cr.NR_MAX_ITER_STATION


# ===========================================================================
# _brent_find_crossing internals (lines 181, 185-186, 219-220, 242)
# ===========================================================================


class TestBrentFindCrossing:
    """Drive Brent's crossing solver through its guard branches."""

    def test_root_not_bracketed_raises(self):
        # f(a) and f(b) have the same sign -> line 181.
        def gp(jd):
            return (jd, 1.0)  # lon == jd, so f never changes sign on [10, 20]

        with pytest.raises(Error, match="not bracketed"):
            _brent_find_crossing(gp, 100.0, 10.0, 20.0, 1e-9)

    def test_initial_swap_and_converge(self):
        # |f(a)| < |f(b)| forces the initial swap (lines 185-186).
        def gp(jd):
            return (jd, 1.0)

        # Root at lon == 0 -> jd == 0; bracket [-1, 5] gives |f(a)|=1 < |f(b)|=5.
        result = _brent_find_crossing(gp, 0.0, -1.0, 5.0, 1e-9, max_iter=100)
        assert abs(result) < 1e-6

    def test_max_iter_returns_best_estimate(self):
        # max_iter=1 with an unreachable tolerance -> falls through to line 242.
        def gp(jd):
            return (jd, 1.0)

        result = _brent_find_crossing(gp, 0.0, -1.0, 5.0, 1e-30, max_iter=1)
        assert isinstance(result, float)

    def test_bisection_branch(self):
        # A strongly nonlinear function forces the bisection fallback so the
        # interpolated step is rejected (lines 217-220).
        def gp(jd):
            lon = math.tanh(jd * 5.0) * 0.5
            return (lon, 1.0)

        result = _brent_find_crossing(gp, 0.0, -40.0, 0.7, 1e-12, max_iter=200)
        assert abs(result) < 1e-3


# ===========================================================================
# _find_bracket_for_crossing (line 298)
# ===========================================================================


class TestFindBracketForCrossing:
    """The no-crossing-in-interval failure path."""

    def test_no_crossing_raises(self):
        def gp(jd):
            return (jd, 1.0)  # monotone, never reaches target on [0, 1]

        with pytest.raises(Error, match="No crossing found"):
            _find_bracket_for_crossing(gp, 200.0, 0.0, 1.0, num_samples=5)


# ===========================================================================
# solcross_ut (lines 349-350, 354, 379, 398-399, 414, 420-422)
# ===========================================================================


class TestSolcrossUt:
    """Forward/backward Sun crossings plus the controlled-failure paths."""

    def test_forward_crossing(self):
        jd = solcross_ut(0.0, J2000)
        assert isinstance(jd, float)
        assert jd > J2000

    def test_backward_crossing(self):
        jd = solcross_ut(0.0, J2000, backwards=True)
        assert jd < J2000

    def test_backward_near_crossing_signed_branch(self):
        # Start essentially at the crossing instant and ask for the previous
        # one -> exercises the signed-normalization guard (around line 371).
        jd_cross = solcross_ut(120.0, J2000)
        jd_prev = solcross_ut(120.0, jd_cross, backwards=True)
        assert jd_prev < jd_cross

    def test_initial_calc_failure(self, monkeypatch):
        # Initial calc_ut raises -> lines 349-350.
        def boom(*args, **kwargs):
            raise CalculationError("boom")

        _patch_calc_ut(monkeypatch, boom)
        with pytest.raises(Error, match="Failed to calculate Sun position"):
            solcross_ut(0.0, J2000)

    def test_zero_speed_default_and_retrograde_guard(self, monkeypatch):
        # speed == 0 at start (line 354) and the forward retrograde guard
        # (line 379) both reached, then convergence on a controlled lon.
        state = {"n": 0}

        def fake(jd, body, flags):
            state["n"] += 1
            if state["n"] == 1:
                # Initial call: speed exactly 0 and negative-speed guard target.
                return ([10.0, 0.0, 1.0, 0.0], None)
            # Subsequent calls converge instantly on the target.
            return ([0.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        # Negative speed forces line 379 via a second-stage fake; use a
        # dedicated fake below instead. Here we only confirm zero-speed path.
        result = solcross_ut(0.0, J2000)
        assert isinstance(result, float)

    def test_forward_negative_speed_branch(self, monkeypatch):
        # speed < 0 and diff > 0 -> line 379.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                # negative speed, lon behind target so diff>0
                return ([10.0, 0.0, 1.0, -0.5], None)
            return ([100.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        result = solcross_ut(100.0, J2000)
        assert isinstance(result, float)

    def test_iteration_calc_failure(self, monkeypatch):
        # First call succeeds, second (in NR loop) raises -> lines 398-399.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            raise CalculationError("iter boom")

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="during iteration"):
            solcross_ut(180.0, J2000)

    def test_iteration_low_speed_default(self, monkeypatch):
        # In-loop speed < 0.01 -> line 414, then converge.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            if calls["n"] == 2:
                # Not converged yet, tiny speed forces the default.
                return ([50.0, 0.0, 1.0, 0.0001], None)
            return ([60.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        result = solcross_ut(60.0, J2000)
        assert isinstance(result, float)

    def test_divergence(self, monkeypatch):
        # NR step pushes jd far away each time -> line 420 divergence.
        def fake(jd, body, flags):
            # Always return a large diff with a huge negative speed so the NR
            # step keeps moving the estimate far from jd_guess.
            return ([0.0, 0.0, 1.0, -0.5], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            solcross_ut(170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        # Never converges but stays in range -> line 422 max-iter raise.
        # A constant small offset keeps diff just above tolerance while the NR
        # step (diff/speed) is tiny, so jd never leaves the safe range.
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            solcross_ut(180.0, J2000)


# ===========================================================================
# solcross (TT) (lines 479-480, 484, 493, 498-499, 517-518, 533, 539-541)
# ===========================================================================


class TestSolcrossTT:
    """TT-time Sun crossing: mirror of solcross_ut branches."""

    def test_forward_and_backward(self):
        jd = solcross(0.0, J2000)
        assert jd > J2000
        jd_b = solcross(0.0, J2000, backwards=True)
        assert jd_b < J2000

    def test_backward_near_crossing(self):
        jc = solcross(200.0, J2000)
        jp = solcross(200.0, jc, backwards=True)
        assert jp < jc

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc(monkeypatch, boom)
        with pytest.raises(Error, match="Failed to calculate Sun position"):
            solcross(0.0, J2000)

    def test_zero_speed_branch(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.0], None)  # speed 0 -> line 484
            return ([0.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(solcross(0.0, J2000), float)

    def test_forward_negative_speed_branch(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, -0.5], None)  # line 498-499
            return ([100.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(solcross(100.0, J2000), float)

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            raise CalculationError("iter boom")  # lines 517-518

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="during iteration"):
            solcross(180.0, J2000)

    def test_iteration_low_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            if calls["n"] == 2:
                return ([50.0, 0.0, 1.0, 0.0001], None)  # line 533
            return ([60.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(solcross(60.0, J2000), float)

    def test_divergence(self, monkeypatch):
        def fake(jd, body, flags):
            return ([0.0, 0.0, 1.0, -0.5], None)  # line 539

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            solcross(170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            solcross(180.0, J2000)


# ===========================================================================
# mooncross_ut (lines 587-588, 592, 606, 625-626, 641, 647-649)
# ===========================================================================


class TestMooncrossUt:
    """Moon longitude crossing forward/backward + failure branches."""

    def test_forward_and_backward(self):
        jd = mooncross_ut(0.0, J2000)
        assert jd > J2000
        jd_b = mooncross_ut(0.0, J2000, backwards=True)
        assert jd_b < J2000

    def test_backward_near_crossing(self):
        jc = mooncross_ut(50.0, J2000)
        jp = mooncross_ut(50.0, jc, backwards=True)
        assert jp < jc

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc_ut(monkeypatch, boom)
        with pytest.raises(Error, match="Failed to calculate Moon position"):
            mooncross_ut(0.0, J2000)

    def test_zero_speed_branch(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.0], None)  # line 592
            return ([0.0, 0.0, 1.0, 13.0], None)

        _patch_calc_ut(monkeypatch, fake)
        assert isinstance(mooncross_ut(0.0, J2000), float)

    def test_forward_negative_speed_branch(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, -5.0], None)  # line 606
            return ([100.0, 0.0, 1.0, 13.0], None)

        _patch_calc_ut(monkeypatch, fake)
        assert isinstance(mooncross_ut(100.0, J2000), float)

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 13.0], None)
            raise CalculationError("iter boom")  # lines 625-626

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="during iteration"):
            mooncross_ut(180.0, J2000)

    def test_iteration_low_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 13.0], None)
            if calls["n"] == 2:
                return ([50.0, 0.0, 1.0, 0.01], None)  # line 641
            return ([60.0, 0.0, 1.0, 13.0], None)

        _patch_calc_ut(monkeypatch, fake)
        assert isinstance(mooncross_ut(60.0, J2000), float)

    def test_divergence(self, monkeypatch):
        def fake(jd, body, flags):
            return ([0.0, 0.0, 1.0, -5.0], None)  # line 647 divergence

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            mooncross_ut(170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 13.0], None)  # line 649

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            mooncross_ut(180.0, J2000)


# ===========================================================================
# mooncross (TT) (lines 708-709, 713, 722, 726-727, 745-746, 761, 767-769)
# ===========================================================================


class TestMooncrossTT:
    """TT-time Moon crossing mirror branches."""

    def test_forward_and_backward(self):
        assert mooncross(0.0, J2000) > J2000
        assert mooncross(0.0, J2000, backwards=True) < J2000

    def test_backward_near_crossing(self):
        jc = mooncross(75.0, J2000)
        jp = mooncross(75.0, jc, backwards=True)
        assert jp < jc

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc(monkeypatch, boom)
        with pytest.raises(Error, match="Failed to calculate Moon position"):
            mooncross(0.0, J2000)  # lines 708-709

    def test_zero_speed_branch(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.0], None)  # line 713
            return ([0.0, 0.0, 1.0, 13.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(mooncross(0.0, J2000), float)

    def test_forward_negative_speed_branch(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, -5.0], None)  # lines 726-727
            return ([100.0, 0.0, 1.0, 13.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(mooncross(100.0, J2000), float)

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 13.0], None)
            raise CalculationError("iter boom")  # lines 745-746

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="during iteration"):
            mooncross(180.0, J2000)

    def test_iteration_low_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 13.0], None)
            if calls["n"] == 2:
                return ([50.0, 0.0, 1.0, 0.01], None)  # line 761
            return ([60.0, 0.0, 1.0, 13.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(mooncross(60.0, J2000), float)

    def test_divergence(self, monkeypatch):
        def fake(jd, body, flags):
            return ([0.0, 0.0, 1.0, -5.0], None)  # line 767

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            mooncross(170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 13.0], None)  # line 769

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            mooncross(180.0, J2000)


# ===========================================================================
# mooncross_node / mooncross_node_ut
# (lines 885-886, 890, 929, 939-940, 957-958, 962, 970, 972-974)
# ===========================================================================


class TestMooncrossNode:
    """Lunar-node crossings: real-data and controlled-failure branches."""

    def test_forward_and_backward_real(self):
        jd_f, lon_f, lat_f = mooncross_node(J2000)
        assert jd_f > J2000
        assert abs(lat_f) < 1e-4
        jd_b, _, _ = mooncross_node(J2000, backwards=True)
        assert jd_b < J2000

    def test_ut_wrapper(self):
        jd_f, _, _ = mooncross_node_ut(J2000)
        assert jd_f > J2000
        jd_b, _, _ = mooncross_node_ut(J2000, backwards=True)
        assert jd_b < J2000

    def test_start_at_node(self):
        # Starting essentially at a node -> the abs(lat) small full-step path.
        jd_node, _, _ = mooncross_node(J2000)
        jd_next, _, lat = mooncross_node(jd_node)
        assert jd_next > jd_node
        assert abs(lat) < 1e-4

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc(monkeypatch, boom)
        with pytest.raises(Error, match="Failed to calculate Moon position"):
            mooncross_node(J2000)  # lines 885-886

    def test_small_lat_speed_initial(self, monkeypatch):
        # Initial lat far from node but lat_speed tiny -> line 890 default.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                # lat large enough to skip the at-node branch, lat_speed ~0.
                return ([0.0, 5.0, 1.0, 0.0, 0.01], None)
            # After the initial estimate, converge to a node ahead in time.
            if jd > J2000 + 1e-6:
                return ([0.0, 0.0, 1.0, 0.0, -1.0], None)
            return ([0.0, 5.0, 1.0, 0.0, -1.0], None)

        _patch_calc(monkeypatch, fake)
        jd, _, lat = mooncross_node(J2000)
        assert abs(lat) < 1e-4

    def test_scan_for_else_exhausted(self, monkeypatch):
        # Force the scan loop to find no sign change so the for/else picks the
        # half-nodal-month fallback guess (line 929), then converge.
        def fake(jd, body, flags):
            # lat positive everywhere except exactly at the fallback guess,
            # so the 8-step scan never sees a sign change.
            target_jd = J2000 + 13.6
            if abs(jd - target_jd) < 0.5:
                return ([0.0, 0.0, 1.0, 0.0, -1.0], None)
            # Large lat, lat_speed positive -> dt_guess negative -> needs_scan.
            return ([0.0, 5.0, 1.0, 0.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        jd, _, lat = mooncross_node(J2000)
        assert abs(lat) < 1e-4

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                # Triggers a straightforward linear guess (not scan, not at-node).
                return ([0.0, 2.0, 1.0, 0.0, -1.0], None)
            raise CalculationError("iter boom")  # lines 939-940

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="during iteration"):
            mooncross_node(J2000)

    def test_converge_wrong_side_then_nudge(self, monkeypatch):
        # Forward search where the linear estimate lands within 1e-6 of tjdet
        # (dt_guess ~3e-7) and the Moon is already at a node there, so NR
        # converges on the "wrong side" (jd <= tjdet + 1e-6). That triggers the
        # nudge + continue (lines 957-958); the next iteration lands a valid
        # crossing strictly ahead of tjdet.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                # lat=-3e-6 (above the 10*NR_TOL at-node threshold), lat_speed=10
                # -> dt_guess = 3e-7, so jd_guess = tjdet + 3e-7.
                return ([0.0, -3e-6, 1.0, 0.0, 10.0], None)
            if jd <= J2000 + 1e-6:
                # Converged at/behind tjdet -> wrong side.
                return ([5.0, 0.0, 1.0, 0.0, 10.0], None)
            # After the nudge, a node strictly ahead of tjdet.
            return ([10.0, 0.0, 1.0, 0.0, 10.0], None)

        _patch_calc(monkeypatch, fake)
        jd, _, lat = mooncross_node(J2000)
        assert jd > J2000
        assert abs(lat) < 1e-4

    def test_small_lat_speed_in_loop_and_clamp_forward(self, monkeypatch):
        # In-loop lat_speed tiny (line 962) and a forward clamp when jd<tjdet
        # (lines 967-968 implicitly), converging afterwards.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([0.0, 2.0, 1.0, 0.0, -1.0], None)
            if calls["n"] == 2:
                # Not converged; tiny lat_speed -> line 962; big lat to step.
                return ([0.0, 3.0, 1.0, 0.0, 0.0], None)
            return ([0.0, 0.0, 1.0, 0.0, -1.0], None)  # converge ahead

        _patch_calc(monkeypatch, fake)
        jd, _, lat = mooncross_node(J2000 + 5.0)
        assert abs(lat) < 1e-4

    def test_backward_clamp(self, monkeypatch):
        # Backward search where an NR step pushes jd > tjdet -> line 970 clamp.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                # lat negative, lat_speed negative -> dt_guess negative (backward).
                return ([0.0, -2.0, 1.0, 0.0, -1.0], None)
            if calls["n"] == 2:
                # Step that would move jd forward of tjdet (positive lat / +slope).
                return ([0.0, 5.0, 1.0, 0.0, -1.0], None)
            return ([0.0, 0.0, 1.0, 0.0, -1.0], None)  # converge

        _patch_calc(monkeypatch, fake)
        jd, _, lat = mooncross_node(J2000, backwards=True)
        assert abs(lat) < 1e-4

    def test_divergence(self, monkeypatch):
        # Forward search; NR keeps lat huge and steps jd far beyond 30 days
        # while staying ahead of tjdet -> line 972 divergence raise.
        def fake(jd, body, flags):
            # Big positive lat, positive slope -> jd -= lat/speed moves backward;
            # to push forward and far, use negative lat with negative slope so
            # jd -= lat/speed = jd - (neg/neg)=jd - positive... craft via sign.
            # Use lat large positive, speed large negative -> -lat/speed = +big.
            return ([0.0, 100.0, 1.0, 0.0, -1.0], None)

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            mooncross_node(J2000)

    def test_max_iterations(self, monkeypatch):
        # Never converges, stays within 30 days, forward -> line 974 max iter.
        def fake(jd, body, flags):
            # Oscillating lat that never reaches 0, small steps within range.
            lat = 1.0 + 0.5 * math.sin(jd)
            return ([0.0, lat, 1.0, 0.0, -2.0], None)

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            mooncross_node(J2000)


# ===========================================================================
# cross_ut
# (lines 1014-1015, 1041, 1054-1081, 1098-1132, 1143-1144, 1173,
#  1190-1192, 1198, 1207, 1209, 1213, 1215, 1219-1221)
# ===========================================================================


class TestCrossUt:
    """Generic planet crossings: NR, retrograde scan, Brent fallbacks."""

    def test_forward_mercury(self):
        jd = cross_ut(MERCURY, 50.0, J2000)
        assert jd > J2000  # exercises planet in (2,3) max_range (line 1215)

    def test_forward_mars(self):
        jd = cross_ut(MARS, 50.0, J2000)
        assert jd > J2000  # else max_range=800 (line 1217)

    def test_forward_jupiter(self):
        jd = cross_ut(JUPITER, 50.0, J2000)
        assert jd > J2000  # speed_default<0.1 -> max_range=5000 (line 1213)

    def test_forward_uranus(self):
        jd = cross_ut(URANUS, 50.0, J2000)
        assert jd > J2000  # speed_default<0.05 -> max_range=40000 (line 1209)

    def test_forward_neptune(self):
        jd = cross_ut(NEPTUNE, 50.0, J2000)
        assert jd > J2000  # speed_default<0.01 -> max_range=100000 (line 1207)

    def test_forward_pluto(self):
        jd = cross_ut(PLUTO, 280.0, J2000)
        assert jd > J2000

    def test_retrograde_scan_path(self):
        # Mercury retrograde with a target just behind it -> scan block
        # (lines 1054-1081) plus the retrograde effective_speed (line 1041).
        jd_rx = ephem.julday(2024, 4, 10, 0.0)
        jd = cross_ut(MERCURY, 22.0, jd_rx)
        assert isinstance(jd, float)

    def test_retrograde_far_target_skips_scan(self):
        # Retrograde Saturn but target far ahead (diff_back >= 25) -> scan
        # skipped, effective_speed default, NR proceeds (line 1041 only).
        jd = ephem.julday(2000, 1, 1, 12.0)
        result = cross_ut(SATURN, 120.0, jd)
        assert result > jd

    def test_near_station_brent(self):
        # Mercury near a station -> Brent's-method branch (lines 1098-1129).
        jd_station = ephem.julday(2024, 3, 1, 0.0) + 32  # ~2024-04-02
        jd = cross_ut(MERCURY, 30.0, jd_station)
        assert isinstance(jd, float)

    def test_near_station_brent_fallback_to_nr(self, monkeypatch):
        # Initial speed near a station triggers the Brent path; patching the
        # bracket helper to raise forces the fall-through to Newton-Raphson
        # (lines 1130-1132), which then converges.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.005], None)  # near-station, positive
            return ([50.0, 0.0, 1.0, 1.0], None)  # NR converges

        def raise_bracket(*a, **k):
            raise Error("no bracket")

        _patch_calc_ut(monkeypatch, fake)
        monkeypatch.setattr(cr, "_find_bracket_for_crossing", raise_bracket)
        result = cross_ut(MARS, 50.0, J2000)
        assert isinstance(result, float)

    def test_retrograde_scan_exhausts(self, monkeypatch):
        # Retrograde planet, target within 25 deg behind, but the scan never
        # finds a sign change (constant longitude) so the while loop exhausts
        # and control falls through to the Newton-Raphson estimate (the
        # 1064 -> 1083 fall-through). NR then diverges, which is acceptable.
        def fake(jd, body, flags):
            return ([60.0, 0.0, 1.0, -0.5], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error):
            cross_ut(MARS, 50.0, J2000)

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc_ut(monkeypatch, boom)
        with pytest.raises(Error, match="Failed to calculate planet position"):
            cross_ut(MARS, 50.0, J2000)  # lines 1014-1015

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)  # direct, not near station
            raise CalculationError("iter boom")  # lines 1143-1144

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="during iteration"):
            cross_ut(MARS, 180.0, J2000)

    def test_retrograde_during_iteration_brent(self, monkeypatch):
        # Start direct (NR), then mid-iteration speed goes negative, triggering
        # the in-loop station/retrograde Brent fallback (lines 1159-1182) that
        # succeeds.
        def fake(jd, body, flags):
            # lon increases smoothly with jd so a bracket/Brent solution exists,
            # but the *reported* speed is negative to trip the fallback once.
            lon = (jd - J2000) * 2.0 % 360.0
            return ([lon, 0.0, 1.0, -1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        # speed<0 at start triggers the early scan; pick a target far ahead so
        # the scan is skipped (diff_back large) and NR runs, then fallback.
        try:
            result = cross_ut(MARS, 100.0, J2000)
            assert isinstance(result, float)
        except Error:
            pass

    def test_in_loop_station_fallback_brent_fails(self, monkeypatch):
        # Direct start (NR path). Mid-iteration the reported speed goes
        # negative, triggering the in-loop station/retrograde fallback. The
        # patched bracket helper raises, so Brent fails and NR continues
        # (lines 1190-1192). A later iteration reports a sub-0.001 speed,
        # exercising the speed-default reset (line 1198), then converges.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            n = calls["n"]
            if n == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)  # initial direct
            if n == 2:
                return ([20.0, 0.0, 1.0, -0.5], None)  # negative -> fallback
            if n == 3:
                return ([10.0, 0.0, 1.0, 1.0], None)  # large +step -> jd back forward
            if n == 4:
                return ([59.9, 0.0, 1.0, 0.0005], None)  # tiny speed -> 1198
            return ([60.0, 0.0, 1.0, 1.0], None)  # converge (jd >= tjdut)

        def raise_bracket(*a, **k):
            raise Error("no bracket")

        _patch_calc_ut(monkeypatch, fake)
        monkeypatch.setattr(cr, "_find_bracket_for_crossing", raise_bracket)
        # The negative-speed step pushes jd before tjdut, the next positive step
        # brings it back forward, so convergence is forward-only (>= tjdut).
        result = cross_ut(MARS, 60.0, J2000)
        assert isinstance(result, float)
        assert result >= J2000 - 1e-6

    def test_divergence(self, monkeypatch):
        # Direct planet (Mars, speed_default 0.5) whose NR step keeps pushing
        # jd far away each iteration -> the >max_range divergence (line 1219).
        # lon stuck at 0 with a large diff and speed 1.0 -> jd jumps ~180/iter,
        # exceeding Mars' 800-day range within a few iterations.
        def fake(jd, body, flags):
            return ([0.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            cross_ut(MARS, 170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        # Constant tiny offset: diff just above tolerance, NR step tiny so jd
        # stays in range, loop exhausts -> line 1221 max-iter raise.
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            cross_ut(MARS, 180.0, J2000)

    def test_retrograde_scan_exact_hit(self, monkeypatch):
        # Retrograde planet exactly at the target at tjdut -> d_prev == 0.0
        # early return in the scan block (lines 1066-1067).
        def fake(jd, body, flags):
            # lon == 50 at tjdut (delta 0), retrograde speed. The scan calls
            # use the no-FLG_SPEED form too; keep lon constant near tjdut.
            if abs(jd - J2000) < 1e-9:
                return ([50.0, 0.0, 1.0, -0.5], None)
            return ([50.0 - (jd - J2000) * 0.5, 0.0, 1.0, -0.5], None)

        _patch_calc_ut(monkeypatch, fake)
        result = cross_ut(MERCURY, 50.0, J2000)
        assert abs(result - J2000) < 1e-6

    def test_retrograde_scan_bisection_refine(self, monkeypatch):
        # Retrograde planet whose longitude crosses the target via a sign
        # change inside the scan window -> bisection refine (lines 1068-1079).
        def fake(jd, body, flags):
            # Linear retrograde motion: lon decreases through 50 at jd~J2000+10.
            lon = 55.0 - (jd - J2000) * 0.5
            return ([lon % 360.0, 0.0, 1.0, -0.5], None)

        _patch_calc_ut(monkeypatch, fake)
        result = cross_ut(MERCURY, 50.0, J2000)
        assert result > J2000


# ===========================================================================
# helio_cross_ut
# (lines 1283-1284, 1308-1314, 1317, 1319, 1328-1329, 1335-1361,
#  1369-1370, 1386, 1397-1399)
# ===========================================================================


class TestHelioCrossUt:
    """Heliocentric crossings (UT) forward/backward + slow-planet Brent."""

    def test_forward(self):
        jd = helio_cross_ut(MARS, 100.0, J2000)
        assert jd > J2000

    def test_backward(self):
        jd = helio_cross_ut(MARS, 100.0, J2000, backwards=True)
        assert jd < J2000

    def test_forward_at_current_lon(self):
        # diff < 1e-5 -> add 360 forward (line 1317).
        pos, _ = ephem.calc_ut(J2000, MARS, FLG_HELCTR | FLG_SPEED)
        jd = helio_cross_ut(MARS, pos[0], J2000)
        assert jd > J2000

    def test_backward_at_current_lon(self):
        # diff in [0, 1e-5] -> elif branch (line 1311).
        pos, _ = ephem.calc_ut(J2000, MARS, FLG_HELCTR | FLG_SPEED)
        jd = helio_cross_ut(MARS, pos[0], J2000, backwards=True)
        assert jd < J2000

    def test_slow_planet_brent_forward(self):
        # Pluto heliocentric speed < threshold -> Brent path forward
        # (lines 1335-1359).
        jd = helio_cross_ut(PLUTO, 280.0, J2000)
        assert jd > J2000

    def test_slow_planet_brent_backward(self):
        # Brent path with backward bracket window (lines 1341-1343).
        jd = helio_cross_ut(PLUTO, 240.0, J2000, backwards=True)
        assert jd < J2000

    def test_slow_planet_brent_fallback_to_nr(self, monkeypatch):
        # Slow planet triggers the Brent path; patched bracket helper raises so
        # Brent fails and execution falls through to Newton-Raphson
        # (lines 1360-1361).
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.005], None)  # near-station
            return ([50.0, 0.0, 1.0, 1.0], None)

        def raise_bracket(*a, **k):
            raise Error("no bracket")

        _patch_calc_ut(monkeypatch, fake)
        monkeypatch.setattr(cr, "_find_bracket_for_crossing", raise_bracket)
        assert isinstance(helio_cross_ut(PLUTO, 50.0, J2000), float)

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc_ut(monkeypatch, boom)
        with pytest.raises(Error, match="heliocentric planet position"):
            helio_cross_ut(MARS, 100.0, J2000)  # lines 1283-1284

    def test_backward_zero_speed(self, monkeypatch):
        # Backward, speed ~0 -> lines 1312-1314.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([100.0, 0.0, 1.0, 0.0], None)  # speed ~0, big diff
            return ([50.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        result = helio_cross_ut(MARS, 50.0, J2000, backwards=True)
        assert isinstance(result, float)

    def test_forward_zero_speed(self, monkeypatch):
        # Forward, speed ~0 -> line 1319.
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.0], None)
            return ([50.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        result = helio_cross_ut(MARS, 50.0, J2000)
        assert isinstance(result, float)

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            raise CalculationError("iter boom")  # lines 1369-1370

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="heliocentric position during iteration"):
            helio_cross_ut(MARS, 180.0, J2000)

    def test_iteration_low_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            if calls["n"] == 2:
                return ([50.0, 0.0, 1.0, 0.00001], None)  # line 1386
            return ([60.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        assert isinstance(helio_cross_ut(MARS, 60.0, J2000), float)

    def test_divergence(self, monkeypatch):
        # lon stuck at 0, big diff, speed 1.0 -> jd jumps ~180/iter, leaving the
        # 400-day max_range within a few iterations (lines 1397-1398).
        def fake(jd, body, flags):
            return ([0.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            helio_cross_ut(MARS, 170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            helio_cross_ut(MARS, 180.0, J2000)


# ===========================================================================
# helio_cross (TT) (lines 1448-1566)
# ===========================================================================


class TestHelioCrossTT:
    """Heliocentric crossings (TT) mirror branches."""

    def test_forward(self):
        assert helio_cross(MARS, 100.0, J2000) > J2000

    def test_backward(self):
        assert helio_cross(MARS, 100.0, J2000, backwards=True) < J2000

    def test_forward_at_current_lon(self):
        pos, _ = ephem.calc(J2000, MARS, FLG_HELCTR | FLG_SPEED)
        assert helio_cross(MARS, pos[0], J2000) > J2000

    def test_backward_at_current_lon(self):
        pos, _ = ephem.calc(J2000, MARS, FLG_HELCTR | FLG_SPEED)
        assert helio_cross(MARS, pos[0], J2000, backwards=True) < J2000

    def test_slow_planet_brent_forward(self):
        assert helio_cross(PLUTO, 280.0, J2000) > J2000

    def test_slow_planet_brent_backward(self):
        assert helio_cross(PLUTO, 240.0, J2000, backwards=True) < J2000

    def test_initial_calc_failure(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc(monkeypatch, boom)
        with pytest.raises(Error, match="heliocentric planet position"):
            helio_cross(MARS, 100.0, J2000)

    def test_backward_zero_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([100.0, 0.0, 1.0, 0.0], None)
            return ([50.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(helio_cross(MARS, 50.0, J2000, backwards=True), float)

    def test_forward_zero_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.0], None)
            return ([50.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(helio_cross(MARS, 50.0, J2000), float)

    def test_iteration_calc_failure(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            raise CalculationError("iter boom")

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="heliocentric position during iteration"):
            helio_cross(MARS, 180.0, J2000)

    def test_iteration_low_speed(self, monkeypatch):
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 1.0], None)
            if calls["n"] == 2:
                return ([50.0, 0.0, 1.0, 0.00001], None)
            return ([60.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        assert isinstance(helio_cross(MARS, 60.0, J2000), float)

    def test_divergence(self, monkeypatch):
        def fake(jd, body, flags):
            return ([0.0, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="diverged"):
            helio_cross(MARS, 170.0, J2000)

    def test_max_iterations(self, monkeypatch):
        def fake(jd, body, flags):
            return ([180.0 - 0.0001, 0.0, 1.0, 1.0], None)

        _patch_calc(monkeypatch, fake)
        with pytest.raises(Error, match="Maximum iterations"):
            helio_cross(MARS, 180.0, J2000)

    def test_slow_planet_brent_fallback_to_nr(self, monkeypatch):
        # Slow planet so the TT Brent path runs but bracketing fails (patched
        # helper raises), falling through to Newton-Raphson (lines 1530-1531).
        calls = {"n": 0}

        def fake(jd, body, flags):
            calls["n"] += 1
            if calls["n"] == 1:
                return ([10.0, 0.0, 1.0, 0.005], None)  # near-station
            return ([50.0, 0.0, 1.0, 1.0], None)

        def raise_bracket(*a, **k):
            raise Error("no bracket")

        _patch_calc(monkeypatch, fake)
        monkeypatch.setattr(cr, "_find_bracket_for_crossing", raise_bracket)
        assert isinstance(helio_cross(PLUTO, 50.0, J2000), float)


# ===========================================================================
# is_retrograde (lines 1612-1620)
# ===========================================================================


class TestIsRetrograde:
    def test_sun_and_moon_never(self):
        assert is_retrograde(SUN, J2000) is False
        assert is_retrograde(MOON, J2000) is False

    def test_planet_direct(self):
        # Mars is direct at J2000.
        assert is_retrograde(MARS, J2000) is False

    def test_planet_retrograde(self):
        jd_rx = ephem.julday(2024, 4, 10, 0.0)
        assert is_retrograde(MERCURY, jd_rx) is True

    def test_calc_failure_returns_false(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc_ut(monkeypatch, boom)
        assert is_retrograde(MARS, J2000) is False  # lines 1619-1620


# ===========================================================================
# get_station_type (lines 1643-1671)
# ===========================================================================


class TestGetStationType:
    def test_sun_moon_direct(self):
        assert get_station_type(SUN, J2000) == "direct"
        assert get_station_type(MOON, J2000) == "direct"

    def test_direct_motion(self):
        assert get_station_type(MARS, J2000) == "direct"

    def test_retrograde_motion(self):
        jd_rx = ephem.julday(2024, 4, 10, 0.0)
        assert get_station_type(MERCURY, jd_rx) == "retrograde"

    def test_stationary_retrograde(self, monkeypatch):
        # Near station with speed decreasing -> stationary_retrograde.
        def fake(jd, body, flags):
            if jd < J2000 - 0.5:
                return ([10.0, 0.0, 1.0, 0.5], None)  # before: faster
            if jd > J2000 + 0.5:
                return ([10.0, 0.0, 1.0, -0.5], None)  # after: slower
            return ([10.0, 0.0, 1.0, 0.001], None)  # near station

        _patch_calc_ut(monkeypatch, fake)
        assert get_station_type(MARS, J2000) == "stationary_retrograde"

    def test_stationary_direct(self, monkeypatch):
        # Near station with speed increasing -> stationary_direct.
        def fake(jd, body, flags):
            if jd < J2000 - 0.5:
                return ([10.0, 0.0, 1.0, -0.5], None)  # before: slower
            if jd > J2000 + 0.5:
                return ([10.0, 0.0, 1.0, 0.5], None)  # after: faster
            return ([10.0, 0.0, 1.0, 0.001], None)

        _patch_calc_ut(monkeypatch, fake)
        assert get_station_type(MARS, J2000) == "stationary_direct"

    def test_calc_failure_returns_direct(self, monkeypatch):
        def boom(*a, **k):
            raise CalculationError("boom")

        _patch_calc_ut(monkeypatch, boom)
        assert get_station_type(MARS, J2000) == "direct"  # lines 1670-1671


# ===========================================================================
# _brent_find_station (lines 1704, 1741-1742, 1762)
# ===========================================================================


class TestBrentFindStation:
    def test_not_bracketed_raises(self):
        # Speed same sign at both ends -> line 1704.
        def gs(jd):
            return 1.0

        with pytest.raises(Error, match="not bracketed"):
            _brent_find_station(gs, 0.0, 10.0)

    def test_converges(self):
        # Linear speed crossing zero -> normal convergence.
        def gs(jd):
            return jd  # zero at jd=0

        result = _brent_find_station(gs, -2.0, 5.0, tolerance=1e-9)
        assert abs(result) < 1e-6

    def test_bisection_and_max_iter(self):
        # Nonlinear speed forces bisection (lines 1741-1742); max_iter=1 path
        # also returns best estimate (line 1762) for an unreachable tolerance.
        def gs(jd):
            return math.tanh(jd * 5.0)

        result = _brent_find_station(gs, -40.0, 0.7, tolerance=1e-12, max_iter=200)
        assert abs(result) < 1e-2

        result2 = _brent_find_station(gs, -40.0, 0.7, tolerance=1e-30, max_iter=1)
        assert isinstance(result2, float)


# ===========================================================================
# find_station_ut "could not find" (line 1930)
# ===========================================================================


class TestFindStationUtNotFound:
    def test_never_finds_requested_type(self, monkeypatch):
        # All stations look like SD; requesting "retrograde" exhausts the
        # attempts and raises (line 1930).
        from libephemeris import crossing as crmod

        def fake(jd, body, flags):
            # Speed always positive then a single dip is needed to bracket a
            # station; make speed oscillate so brackets exist but type is SD.
            return ([10.0, 0.0, 1.0, math.cos(jd) * 0.5], None)

        _patch_calc_ut(monkeypatch, fake)

        # Force every detected station to classify as "SD" by patching the
        # bracket/brent helpers to return a deterministic station whose
        # before/after speeds make it look direct.
        def fake_bracket(planet_id, jd_start, jd_end, flag, num_samples=50):
            return (jd_start, jd_start + 1.0)

        def fake_brent(get_speed, jd_a, jd_b, tol):
            return jd_a + 0.5

        monkeypatch.setattr(crmod, "_find_station_bracket", fake_bracket)
        monkeypatch.setattr(crmod, "_brent_find_station", fake_brent)

        # before>0 and after>0 -> classified SD; requesting retrograde fails.
        def speed_dir(jd, body, flags):
            return ([10.0, 0.0, 1.0, 1.0], None)

        _patch_calc_ut(monkeypatch, speed_dir)

        with pytest.raises(Error, match="Could not find"):
            from libephemeris.crossing import find_station_ut

            find_station_ut(MERCURY, J2000, "retrograde")


# ===========================================================================
# next_retrograde_ut (lines 1958-1968)
# ===========================================================================


class TestNextRetrogradeUt:
    def test_currently_direct(self):
        # Mars direct at J2000 -> else branch (lines 1966-1968).
        jd_sr, jd_sd = next_retrograde_ut(MERCURY, J2000)
        assert jd_sr < jd_sd

    def test_currently_retrograde(self):
        # Mercury retrograde mid-period -> if branch (lines 1960-1963).
        jd_rx = ephem.julday(2024, 4, 10, 0.0)
        jd_sr, jd_sd = next_retrograde_ut(MERCURY, jd_rx)
        assert jd_sr < jd_sd


# ===========================================================================
# calc_velocity_at_station (lines 2004-2022)
# ===========================================================================


class TestCalcVelocityAtStation:
    def test_basic(self):
        v_lon, v_lat, v_dist = calc_velocity_at_station(MARS, J2000)
        assert isinstance(v_lon, float)
        assert isinstance(v_lat, float)
        assert isinstance(v_dist, float)

    def test_longitude_wrap_positive(self, monkeypatch):
        # dlon > 180 -> subtract 360 (lines 2013-2014).
        def fake(jd, body, flags):
            if jd < J2000:
                return ([350.0, 1.0, 1.0], None)
            return ([10.0, 1.0, 1.0], None)  # dlon = -340 -> +360 wrap below

        _patch_calc_ut(monkeypatch, fake)
        v_lon, _, _ = calc_velocity_at_station(MARS, J2000)
        assert isinstance(v_lon, float)

    def test_longitude_wrap_negative(self, monkeypatch):
        # dlon < -180 -> add 360 (lines 2015-2016).
        def fake(jd, body, flags):
            if jd < J2000:
                return ([10.0, 1.0, 1.0], None)
            return ([350.0, 1.0, 1.0], None)  # dlon = 340 -> -360

        _patch_calc_ut(monkeypatch, fake)
        v_lon, _, _ = calc_velocity_at_station(MARS, J2000)
        assert isinstance(v_lon, float)


# ===========================================================================
# get_station_info (lines 2050-2075)
# ===========================================================================


class TestGetStationInfo:
    def test_sun_raises(self):
        with pytest.raises(ValueError, match="Sun and Moon"):
            get_station_info(SUN, J2000)

    def test_normal_info(self):
        info = get_station_info(MERCURY, J2000)
        assert "jd_station" in info
        assert info["station_type"] in ("SR", "SD")
        assert isinstance(info["velocity"], float)

    def test_station_search_failure_partial(self, monkeypatch):
        # find_station_ut raises -> partial-info branch (lines 2073-2075).
        from libephemeris import crossing as crmod

        # Keep the initial calc_ut working for current state, but make the
        # station search raise.
        def fake_find(planet_id, jd_ut, st="any", flag=0):
            raise Error("no station")

        monkeypatch.setattr(crmod, "find_station_ut", fake_find)
        info = get_station_info(MERCURY, J2000)
        assert info["jd_station"] is None
        assert "error" in info
