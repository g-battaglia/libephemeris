"""Defining-condition invariants for the event solvers (oracle-free).

Each finder's result must satisfy the condition it solves for — assert it
directly (no external truth needed): a crossing's longitude equals the target, a
station's speed is ~0 and the motion reverses, a node crossing has zero latitude,
and forward/backward searches land on the correct side. This is the property that
caught the cross_ut bugs; here it guards the rest of the family.
"""
from __future__ import annotations

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN, MOON, MERCURY, VENUS, MARS, JUPITER, SATURN, FLG_SWIEPH, FLG_SPEED,
    FLG_HELCTR,
)

_PLANETS = [MERCURY, VENUS, MARS, JUPITER, SATURN]
_STARTS = [2451545.0 + d for d in (-2500.0, -800.0, 0.0, 1200.0, 3500.0)]


def _wrap(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


def _speed(body, jd):
    return L.calc_ut(jd, body, FLG_SWIEPH | FLG_SPEED)[0][3]


@pytest.fixture(autouse=True)
def _skyfield():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")


@pytest.mark.slow
@pytest.mark.parametrize("body", _PLANETS)
def test_find_station_is_zero_speed_reversal(body):
    for start in _STARTS:
        ts, kind = L.find_station_ut(body, start)
        assert ts >= start - 1e-6
        assert abs(_speed(body, ts)) < 0.02, f"{body}: speed {_speed(body, ts)}"
        # motion reverses across the station
        assert (_speed(body, ts - 1.0) < 0) != (_speed(body, ts + 1.0) < 0)
        assert kind in ("SR", "SD")


@pytest.mark.slow
@pytest.mark.parametrize("body", _PLANETS)
def test_next_retrograde_enters_retrograde(body):
    for start in _STARTS:
        rstart, rend = L.next_retrograde_ut(body, start)
        assert rstart >= start - 1e-6 and rend > rstart
        assert abs(_speed(body, rstart)) < 0.03
        assert _speed(body, rstart + 0.5) < 0.01  # retrograde just after start


@pytest.mark.slow
@pytest.mark.parametrize("body", _PLANETS)
def test_helio_cross_forward_and_backward(body):
    for start in _STARTS:
        for target in (37.0, 211.0):
            fwd = L.helio_cross_ut(body, target, start, FLG_SWIEPH)
            assert fwd >= start - 1e-6
            assert _wrap(L.calc_ut(fwd, body, FLG_SWIEPH | FLG_HELCTR)[0][0], target) < 1e-3
            bwd = L.helio_cross_ut(body, target, start, FLG_SWIEPH, True)
            assert bwd <= start + 1e-6
            assert _wrap(L.calc_ut(bwd, body, FLG_SWIEPH | FLG_HELCTR)[0][0], target) < 1e-3


@pytest.mark.slow
def test_mooncross_node_zero_latitude():
    for start in _STARTS:
        jd, _lon, _lat = L.mooncross_node_ut(start, FLG_SWIEPH)
        assert jd >= start - 1e-6
        assert abs(L.calc_ut(jd, MOON, FLG_SWIEPH)[0][1]) * 3600.0 < 1.0  # |lat| < 1"
        bwd = L.mooncross_node_ut(start, FLG_SWIEPH, True)[0]
        assert bwd <= start + 1e-6


@pytest.mark.slow
@pytest.mark.parametrize("fn,body", [(L.solcross_ut, SUN), (L.mooncross_ut, MOON)])
def test_solcross_mooncross_backward(fn, body):
    for start in _STARTS:
        for target in (15.0, 200.0, 333.0):
            bwd = fn(target, start, FLG_SWIEPH, True)
            assert bwd <= start + 1e-6
            assert _wrap(L.calc_ut(bwd, body, FLG_SWIEPH)[0][0], target) < 1e-3
