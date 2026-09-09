# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Analytical light-time solution for uniform rectilinear motion."""

from __future__ import annotations

import math

import pytest

from libephemeris import fast_calc as fc
from libephemeris.constants import EARTH, FLG_NOABERR, FLG_NOGDEFL, SUN


@pytest.mark.parametrize("jd", [2415020.0, 2451545.0, 2490000.125])
@pytest.mark.parametrize("speed", [-0.01, 0.01])
def test_affine_light_time_position_and_velocity(monkeypatch, jd, speed):
    """r(t-lt)=r(t)/(1+v/c) for a stationary observer and radial motion."""
    distance = 0.00257

    class Reader:
        _bodies = {}

        def has_body(self, body_id):
            return body_id in (EARTH, SUN)

        def eval_body(self, body_id, epoch):
            return self._eval_body_split(body_id, epoch, 0.0)

        def _eval_body_split(self, body_id, epoch, offset):
            if body_id == EARTH:
                return (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
            dt = math.fsum((epoch, -jd, offset))
            return (distance + speed * dt, 0.0, 0.0), (speed, 0.0, 0.0)

    # Isolate light time from frame transformations and optical reductions.
    monkeypatch.setattr(fc, "_frame_transform", lambda vector, *_args: vector)
    position = fc._pipeline_icrs(
        Reader(), jd, SUN, FLG_NOABERR | FLG_NOGDEFL, want_xyz=True
    )
    state = fc._pipeline_icrs(
        Reader(),
        jd,
        SUN,
        FLG_NOABERR | FLG_NOGDEFL,
        want_xyz=True,
        want_velocity=True,
    )
    ratio = speed / fc.C_LIGHT_AU_DAY
    expected_position = distance / (1.0 + ratio)
    expected_speed = speed / (1.0 + ratio)
    # Three fixed-point iterations: a known geometric remainder, not a fitted limit.
    position_bound = distance * abs(ratio) ** 4 / (1.0 - abs(ratio))
    position_bound += 8.0 * math.ulp(distance)
    assert abs(position[0] - expected_position) <= position_bound
    assert position[1:] == (0.0, 0.0)
    assert state[:3] == position
    # The published velocity passes through a stencil, so include its roundoff.
    speed_bound = 8.0 * math.ulp(distance) / fc._VEL_H
    assert abs(state[3] - expected_speed) <= speed_bound
    assert state[4:] == (0.0, 0.0)
