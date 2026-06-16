"""Tests for backward crossing search.

Covers ``backwards=True`` in ``solcross_ut``, ``mooncross_ut``,
``mooncross_node_ut`` (and their TT variants), added for planetary
return navigation where users step backward through past returns.
"""

from __future__ import annotations

import pytest
import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    FLG_SWIEPH,
)

JD_J2000 = 2451545.0

# Nominal cycle lengths used for sanity bounds
TROPICAL_YEAR = 365.2422
SIDEREAL_MONTH = 27.3217
NODAL_HALF_MONTH = 13.61


@pytest.mark.unit
class TestSolcrossBackwards:
    """Backward Sun crossings."""

    @pytest.mark.parametrize("degree", [0.0, 90.0, 180.0, 270.0])
    def test_solcross_backwards_within_year(self, degree):
        """Previous Sun crossing must be within the past year."""
        jd = swe.solcross_ut(degree, JD_J2000, FLG_SWIEPH, backwards=True)
        assert jd < JD_J2000
        assert JD_J2000 - jd < 370, f"Previous {degree}° too far back: {JD_J2000 - jd} days"
        sun, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
        err = abs(sun[0] - degree)
        if err > 180:
            err = 360 - err
        assert err < 0.001, f"Sun at {sun[0]}° but expected {degree}°"

    def test_solcross_round_trip_from_crossing(self):
        """Backward search from a crossing yields the previous one ~1 year earlier."""
        jd_next = swe.solcross_ut(55.0, JD_J2000)
        jd_prev = swe.solcross_ut(55.0, jd_next, backwards=True)
        delta = jd_next - jd_prev
        assert abs(delta - TROPICAL_YEAR) < 2.0, f"Round-trip delta {delta} not ~1 yr"

    def test_solcross_forward_back_symmetric(self):
        """next then back yields the previous crossing; fwd from there = next.

        A forward search started exactly AT a crossing returns that
        same instant (pyswisseph-verified: swe.solcross_ut from its own
        returned crossing gives dt = 0.000000000); navigation steps
        forward by adding a small epsilon past the found crossing.
        """
        jd_next = swe.solcross_ut(55.0, JD_J2000)
        jd_prev = swe.solcross_ut(55.0, jd_next, backwards=True)
        # exactly at a crossing: the crossing "now" is returned
        jd_same = swe.solcross_ut(55.0, jd_prev)
        assert abs(jd_same - jd_prev) * 86400 < 0.1
        # epsilon-forward: the next crossing
        jd_forward_again = swe.solcross_ut(55.0, jd_prev + 1.0 / 86400.0)
        err_sec = abs(jd_forward_again - jd_next) * 86400
        assert err_sec < 0.1, f"Forward-from-prev ≠ next: {err_sec}s"

    def test_solcross_multi_step_navigation(self):
        """Multiple forward steps + same number backward return to the first crossing."""
        jd = JD_J2000
        forward_marks = []
        for _ in range(3):
            jd = swe.solcross_ut(55.0, jd)
            forward_marks.append(jd)
            # Epsilon-forward past this crossing so the next search advances:
            # like pyswisseph, solcross_ut from an exact crossing returns that
            # same crossing (see test_solcross_forward_idempotent above).
            jd += 1.0 / 86400.0
        # Navigate backward from the last mark exactly: the backward at-crossing
        # guard steps a full cycle to the previous crossing.
        jd = forward_marks[-1]
        for expected in reversed(forward_marks[:-1]):
            jd = swe.solcross_ut(55.0, jd, backwards=True)
            assert abs(jd - expected) * 86400 < 0.1


@pytest.mark.unit
class TestMooncrossBackwards:
    """Backward Moon crossings."""

    @pytest.mark.parametrize("degree", [0.0, 90.0, 180.0, 270.0])
    def test_mooncross_backwards_within_month(self, degree):
        """Previous Moon crossing must be within the past sidereal month."""
        jd = swe.mooncross_ut(degree, JD_J2000, FLG_SWIEPH, backwards=True)
        assert jd < JD_J2000
        assert JD_J2000 - jd < 30, f"Previous {degree}° too far back: {JD_J2000 - jd} days"
        moon, _ = swe.calc_ut(jd, MOON, FLG_SWIEPH)
        err = abs(moon[0] - degree)
        if err > 180:
            err = 360 - err
        assert err < 0.001, f"Moon at {moon[0]}° but expected {degree}°"

    def test_mooncross_round_trip(self):
        """Next-then-back from a crossing yields a crossing ~1 sidereal month earlier."""
        jd_next = swe.mooncross_ut(120.0, JD_J2000)
        jd_prev = swe.mooncross_ut(120.0, jd_next, backwards=True)
        delta = jd_next - jd_prev
        assert abs(delta - SIDEREAL_MONTH) < 0.5, f"Delta {delta} ≠ ~27.32d"

    def test_mooncross_multi_step(self):
        """Forward/backward navigation with epsilon-stepping past each
        found crossing (a search started exactly at a crossing may
        return that same instant — see the solcross test above)."""
        eps = 1.0 / 86400.0
        jd = JD_J2000
        forward_marks = []
        for _ in range(3):
            jd = swe.mooncross_ut(120.0, jd + eps)
            forward_marks.append(jd)
        for expected in reversed(forward_marks[:-1]):
            jd = swe.mooncross_ut(120.0, jd - eps, backwards=True)
            assert abs(jd - expected) * 86400 < 0.1


@pytest.mark.unit
class TestMooncrossNodeBackwards:
    """Backward Moon-node crossings."""

    def test_mooncross_node_backwards_within_half_month(self):
        """Previous node crossing must be in [JD0 - 14d, JD0)."""
        jd, _, lat = swe.mooncross_node_ut(JD_J2000, FLG_SWIEPH, backwards=True)
        assert jd < JD_J2000
        assert JD_J2000 - jd < NODAL_HALF_MONTH + 2
        assert abs(lat) < 1e-4

    def test_mooncross_node_round_trip(self):
        jd_next, _, _ = swe.mooncross_node_ut(JD_J2000)
        jd_prev, _, _ = swe.mooncross_node_ut(jd_next, backwards=True)
        delta = jd_next - jd_prev
        assert abs(delta - NODAL_HALF_MONTH) < 1.5

    def test_mooncross_node_multi_step(self):
        jd = JD_J2000
        forward_marks = []
        for _ in range(3):
            jd, _, _ = swe.mooncross_node_ut(jd)
            forward_marks.append(jd)
        for expected in reversed(forward_marks[:-1]):
            jd, _, _ = swe.mooncross_node_ut(jd, backwards=True)
            assert abs(jd - expected) * 86400 < 0.1


@pytest.mark.unit
class TestTTVariantsBackwards:
    """TT (Terrestrial Time) variants must also support backwards."""

    def test_solcross_tt_backwards(self):
        jd = swe.solcross(0.0, JD_J2000, FLG_SWIEPH, backwards=True)
        assert jd < JD_J2000

    def test_mooncross_tt_backwards(self):
        jd = swe.mooncross(0.0, JD_J2000, FLG_SWIEPH, backwards=True)
        assert jd < JD_J2000

    def test_mooncross_node_tt_backwards(self):
        jd, _, lat = swe.mooncross_node(JD_J2000, FLG_SWIEPH, backwards=True)
        assert jd < JD_J2000
        assert abs(lat) < 1e-4


@pytest.mark.unit
class TestJustAfterCrossing:
    """Backward search from a point moments after a crossing must return
    that same crossing (seconds ago), not jump a full cycle back.

    Regression test: earlier implementations used a ~1e-3° dead-zone guard
    on solar/lunar and a 0.001-day (86 s) backward-only rejection window
    on the node crossing. Both produced wrong results when the caller was
    within ~minute of the event.
    """

    # Seconds → JD fraction
    _SEC = 1.0 / 86400.0

    @pytest.mark.parametrize("offset_sec", [0.5, 1.0, 10.0, 60.0, 80.0])
    def test_solcross_just_after(self, offset_sec):
        """Sub-second to 80 s after a 55° crossing, previous must be that same crossing."""
        jd_cross = swe.solcross_ut(55.0, JD_J2000)
        jd_start = jd_cross + offset_sec * self._SEC
        jd_prev = swe.solcross_ut(55.0, jd_start, backwards=True)
        # Must land on the same crossing, not the one a full year earlier
        err_sec = abs(jd_prev - jd_cross) * 86400
        assert err_sec < 5.0, (
            f"Previous crossing off by {err_sec:.2f} s from the one {offset_sec}s ago "
            f"— suggests the backward dead-zone jumped a full cycle."
        )

    @pytest.mark.parametrize("offset_sec", [0.5, 1.0, 10.0, 60.0, 80.0])
    def test_mooncross_just_after(self, offset_sec):
        jd_cross = swe.mooncross_ut(120.0, JD_J2000)
        jd_start = jd_cross + offset_sec * self._SEC
        jd_prev = swe.mooncross_ut(120.0, jd_start, backwards=True)
        err_sec = abs(jd_prev - jd_cross) * 86400
        assert err_sec < 5.0, (
            f"Previous moon crossing off by {err_sec:.2f} s from the one {offset_sec}s ago"
        )

    @pytest.mark.parametrize("offset_sec", [0.5, 1.0, 10.0, 60.0, 80.0])
    def test_mooncross_node_just_after(self, offset_sec):
        """Sub-second to 80 s after a node crossing, previous must be that same crossing."""
        jd_cross, _, _ = swe.mooncross_node_ut(JD_J2000)
        jd_start = jd_cross + offset_sec * self._SEC
        jd_prev, _, lat = swe.mooncross_node_ut(jd_start, backwards=True)
        err_sec = abs(jd_prev - jd_cross) * 86400
        assert err_sec < 5.0, (
            f"Previous node crossing off by {err_sec:.2f} s — the 0.001-day "
            f"backward rejection window likely bounced us half a nodal month."
        )
        assert abs(lat) < 1e-4

    def test_solcross_exactly_at_crossing_steps_full_cycle(self):
        """Calling backwards=True at the exact crossing JD must yield the previous cycle."""
        jd_cross = swe.solcross_ut(55.0, JD_J2000)
        jd_prev = swe.solcross_ut(55.0, jd_cross, backwards=True)
        delta_days = jd_cross - jd_prev
        assert abs(delta_days - TROPICAL_YEAR) < 2.0, (
            f"Backward from exact crossing must land ~1 tropical year earlier; got {delta_days} d"
        )

    def test_mooncross_exactly_at_crossing_steps_full_cycle(self):
        jd_cross = swe.mooncross_ut(120.0, JD_J2000)
        jd_prev = swe.mooncross_ut(120.0, jd_cross, backwards=True)
        delta_days = jd_cross - jd_prev
        assert abs(delta_days - SIDEREAL_MONTH) < 0.5

    def test_mooncross_node_exactly_at_crossing_steps_half_nodal_month(self):
        jd_cross, _, _ = swe.mooncross_node_ut(JD_J2000)
        jd_prev, _, _ = swe.mooncross_node_ut(jd_cross, backwards=True)
        delta_days = jd_cross - jd_prev
        assert abs(delta_days - NODAL_HALF_MONTH) < 1.5
