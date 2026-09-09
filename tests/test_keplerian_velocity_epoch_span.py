# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Derivatives use the spacing between the epochs actually evaluated."""

from __future__ import annotations

import math
from types import SimpleNamespace

import pytest

from libephemeris import constants as C
from libephemeris import planets


@pytest.mark.parametrize("epoch", [2396763.5, 2451545.0, -3027092.5, 8000000.0])
def test_linear_motion_uses_actual_epoch_span(monkeypatch, epoch):
    evaluations = []

    def position(jd, body, flags, ephemeris):
        evaluations.append(jd)
        dt = jd - epoch
        return 100.0 + dt, 10.0 - 2.0 * dt, 5.0 + 0.5 * dt

    monkeypatch.setattr(planets, "_keplerian_position_at", position)
    monkeypatch.setattr(
        planets, "_maybe_equatorial_convert", lambda result, *_args: result
    )
    flags = C.FLG_HELCTR | C.FLG_SPEED
    result, returned_flags = planets._calc_keplerian_fallback(
        SimpleNamespace(tt=epoch), C.CHIRON, flags, None
    )

    assert len(evaluations) == 3
    assert evaluations[1] < epoch < evaluations[2]
    assert result[:3] == (100.0, 10.0, 5.0)
    assert result[3:] == (1.0, -2.0, 0.5)
    assert returned_flags == flags


@pytest.mark.parametrize("direction", [-1.0, 1.0])
def test_velocity_wraps_longitude_in_both_directions(monkeypatch, direction):
    epoch = 2451545.0

    def position(jd, *_args):
        return (direction * (jd - epoch)) % 360.0, 0.0, 1.0

    monkeypatch.setattr(planets, "_keplerian_position_at", position)
    monkeypatch.setattr(
        planets, "_maybe_equatorial_convert", lambda result, *_args: result
    )
    result, _flags = planets._calc_keplerian_fallback(
        SimpleNamespace(tt=epoch), C.CHIRON, C.FLG_HELCTR | C.FLG_SPEED, None
    )
    span = (epoch + 1.0 / 86400.0) - (epoch - 1.0 / 86400.0)
    # This bounds roundoff while unwrapping a full turn, not baseline data.
    assert abs(result[3] - direction) <= 2.0 * math.ulp(360.0 / span)
    assert result[:3] == (0.0, 0.0, 1.0)
    assert result[4:] == (0.0, 0.0)


def test_no_speed_evaluates_only_the_central_epoch(monkeypatch):
    epoch = 2451545.0
    calls = []

    def position(jd, *_args):
        calls.append(jd)
        return 20.0, 3.0, 2.0

    monkeypatch.setattr(planets, "_keplerian_position_at", position)
    monkeypatch.setattr(
        planets, "_maybe_equatorial_convert", lambda result, *_args: result
    )
    result, flags = planets._calc_keplerian_fallback(
        SimpleNamespace(tt=epoch), C.CHIRON, C.FLG_HELCTR, None
    )
    assert calls == [epoch]
    assert result == (20.0, 3.0, 2.0, 0.0, 0.0, 0.0)
    assert flags == C.FLG_HELCTR
