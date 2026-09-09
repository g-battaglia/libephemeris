# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Axis distance from exactly representable synthetic states."""

from __future__ import annotations

import math

import pytest

from libephemeris import eclipse
from libephemeris.constants import MOON, SUN


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("offset", [0.0, 2.0**-40, -(2.0**-40), 2.0**-20])
def test_small_perpendicular_distance_is_not_lost(monkeypatch, axis, offset):
    """Distance from a Cartesian axis is the magnitude of the sole component."""
    sun = [0.0, 0.0, 0.0]
    moon = [0.0, 0.0, 0.0]
    sun[axis] = -1.0
    moon[axis] = 2.0**-9
    moon[(axis + 1) % 3] = offset

    def calc_ut(jd, body, flags):
        assert body in (SUN, MOON)
        vector = sun if body == SUN else moon
        return tuple(vector) + (0.0, 0.0, 0.0), flags

    monkeypatch.setattr(eclipse, "calc_ut", calc_ut)
    _word, _attr, geometry = eclipse._lun_how_core(2451545.0)

    assert geometry.axis_offset_km == abs(offset) * eclipse._ECL_AU_KM
    geometry.validate()


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("offset", [0.0, 2.0**-40, -(2.0**-40), 2.0**-20])
def test_small_opposition_angle_is_not_lost(monkeypatch, axis, offset):
    """The synthetic right triangle determines the angle without arccos."""
    sun = [0.0, 0.0, 0.0]
    moon = [0.0, 0.0, 0.0]
    axial_distance = 2.0**-9
    sun[axis] = -1.0
    moon[axis] = axial_distance
    moon[(axis + 1) % 3] = offset

    def calc_ut(jd, body, flags):
        assert body in (SUN, MOON)
        vector = sun if body == SUN else moon
        return tuple(vector) + (0.0, 0.0, 0.0), flags

    monkeypatch.setattr(eclipse, "calc_ut", calc_ut)
    word, attr, _geometry = eclipse._lun_how_core(2451545.0)

    assert word != 0
    expected = math.degrees(math.atan(abs(offset) / axial_distance))
    assert attr[7] == expected
