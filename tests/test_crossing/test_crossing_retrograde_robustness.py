"""Regression tests for cross_ut robustness near retrograde loops.

cross_ut must return the *first* time the planet's ecliptic longitude crosses the
target, searching strictly forward. Three real bugs were found where it did not:

  1. behind-target reachable via an imminent retrograde dip -> returned a crossing
     up to ~60 years late (Saturn: +21773 d where the truth was +40 d);
  2. retrograde-now with a recently-passed (apparently-ahead) target -> Newton
     stepped backward and returned a crossing BEFORE tjdut (forward-only broken);
  3. targets near a station -> the search diverged and raised instead of finding
     the existing crossing.

The previous tests only asserted ``isinstance(result, float)``, so none of these
were caught. These assert the crossing against an independent fine scan of the
library's own longitude curve, plus the forward-only invariant.
"""
from __future__ import annotations

import pytest

import libephemeris as L
from libephemeris.constants import (
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    FLG_SWIEPH,
)


def _lon(planet: int, jd: float) -> float:
    return L.calc_ut(jd, planet, FLG_SWIEPH)[0][0]


def _first_crossing(planet: int, jd_start: float, target: float,
                    max_days: float = 950.0, step: float = 0.5):
    """Independent first genuine sign-change crossing of the library's longitude
    curve (the reference cross_ut should reproduce). None if none in window."""
    prev = _lon(planet, jd_start)
    jd = jd_start
    while jd < jd_start + max_days:
        jd += step
        cur = _lon(planet, jd)
        d1 = ((target - prev + 180.0) % 360.0) - 180.0
        d2 = ((target - cur + 180.0) % 360.0) - 180.0
        if (d1 < 0) != (d2 < 0) and abs(d2 - d1) < 180.0:
            return jd - step / 2.0
        prev = cur
    return None


@pytest.fixture(autouse=True)
def _skyfield_mode():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")


def test_cross_ut_is_forward_only():
    """cross_ut must never return a time before tjdut — including when the planet
    is retrograde and the target lies just 'ahead' (recently passed). Bug #2."""
    jd0 = L.julday(2020, 8, 1, 0.0, L.GREG_CAL)
    for di in (56, 64, 72, 80, 88):  # across Mars' 2020 retrograde
        start = jd0 + di
        base = _lon(MARS, start)
        target = (base + 3.0) % 360.0  # 3 deg "ahead" while retrograde
        jd = L.cross_ut(MARS, target, start, FLG_SWIEPH)
        assert jd >= start - 1e-6, f"cross_ut went backward: {jd - start:.2f} d"
        truth = _first_crossing(MARS, start, target)
        assert truth is not None
        assert abs(jd - truth) < 2.0, f"di={di}: {jd-start:.1f} vs {truth-start:.1f}"


def test_cross_ut_far_behind_prograde_is_forward():
    """A target far behind a prograde planet (beyond retrograde reach) is first
    crossed a WHOLE ORBIT ahead, not in the past. Bug #4 (v7 fuzz): Newton
    converged to / oscillated around the past crossing and returned a negative
    time (or raised); the forward bracket fallback must return the future one.
    Jupiter ~+11 yr, Uranus ~+72 yr ahead.
    """
    for body, start, off in ((JUPITER, 2447651.3, -20.6), (URANUS, 2447883.3, -48.8)):
        base = _lon(body, start)
        target = (base + off) % 360.0
        jd = L.cross_ut(body, target, start, FLG_SWIEPH)
        assert jd >= start - 1e-6, f"body {body} went backward: {jd - start:.0f} d"
        d = abs(((_lon(body, jd) - target + 180.0) % 360.0) - 180.0)
        assert d < 0.01, f"body {body}: wrong crossing longitude {d:.4f} deg"


def test_cross_ut_behind_target_not_decades_late():
    """A target just behind a (still prograde) planet approaching its station is
    reached via the imminent retrograde dip — not a full orbit later. Bug #1."""
    start = L.julday(2020, 5, 1, 0.0, L.GREG_CAL)  # Saturn approaching retro
    target = 301.2596  # ~0.6 deg behind Saturn here
    jd = L.cross_ut(SATURN, target, start, FLG_SWIEPH)
    assert jd - start < 90.0, f"returned {jd-start:.0f} d late (decades-late bug)"
    truth = _first_crossing(SATURN, start, target)
    assert truth is not None and abs(jd - truth) < 2.0


def test_cross_ut_near_station_does_not_diverge():
    """Targets near a retrograde station must be found, not raise 'diverged'.
    Bug #3 (early-September 2020 Mars station)."""
    jd0 = L.julday(2020, 9, 1, 0.0, L.GREG_CAL)
    for di in (1, 4, 7):
        start = jd0 + di
        base = _lon(MARS, start)
        for off in (8.0, 20.0):
            target = (base + off) % 360.0
            jd = L.cross_ut(MARS, target, start, FLG_SWIEPH)  # must not raise
            assert jd >= start - 1e-6


@pytest.mark.slow
@pytest.mark.parametrize("planet,year,month", [
    (MARS, 2020, 8),
    (SATURN, 2020, 3),
])
def test_cross_ut_behind_targets_genuine(planet, year, month):
    """Breadth guard: targets a few degrees BEHIND the planet across its
    retrograde window are reached via the imminent retrograde dip (forward,
    within ~one synodic period) and reproduce an independent fine scan — never a
    full orbit / decades late. (The exhaustive cross-planet sweep with proper
    tangent handling lives in the local validation harness.)"""
    jd0 = L.julday(year, month, 1, 0.0, L.GREG_CAL)
    checked = 0
    for di in range(0, 140, 10):
        start = jd0 + di
        base = _lon(planet, start)
        for off in (-1.0, -3.0, -6.0):  # behind targets -> retrograde crossing
            target = (base + off) % 360.0
            truth = _first_crossing(planet, start, target)
            if truth is None or truth - start > 700.0:
                continue
            jd = L.cross_ut(planet, target, start, FLG_SWIEPH)
            assert jd >= start - 1e-6
            assert abs(jd - truth) < 2.0, (
                f"{planet} di={di} off={off}: {jd-start:.1f} vs {truth-start:.1f}"
            )
            checked += 1
    assert checked > 0
