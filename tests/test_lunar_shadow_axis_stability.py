# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Distanza dall’asse con stati sintetici rappresentabili esattamente."""

from __future__ import annotations

import pytest

from libephemeris import eclipse
from libephemeris.constants import MOON, SUN


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("offset", [0.0, 2.0**-40, -(2.0**-40), 2.0**-20])
def test_small_perpendicular_distance_is_not_lost(monkeypatch, axis, offset):
    """La distanza da un asse cartesiano è il modulo dell’unica componente."""
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
