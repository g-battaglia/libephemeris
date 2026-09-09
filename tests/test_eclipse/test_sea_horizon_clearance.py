# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""La richiesta di orizzonte marino applica la separazione API dichiarata."""

from __future__ import annotations

import pytest

import libephemeris as le
from libephemeris.eclipse import _SEA_HORIZON_CLEARANCE_DEG
from libephemeris.refraction import calc_dip


@pytest.mark.parametrize("height", [0.0, 1000.0])
@pytest.mark.parametrize("event", [le.CALC_RISE, le.CALC_SET])
def test_sea_horizon_matches_explicit_dip_and_clearance(height, event):
    assert _SEA_HORIZON_CLEARANCE_DEG == 0.0001
    geopos = (12.5, 41.9, height)
    horizon = calc_dip(height, atpress=1013.25, attemp=15.0)
    horizon += _SEA_HORIZON_CLEARANCE_DEG
    sentinel = le.rise_trans_true_hor(
        2451545.0, le.SUN, event, geopos, 1013.25, 15.0, -100.0
    )
    explicit = le.rise_trans_true_hor(
        2451545.0, le.SUN, event, geopos, 1013.25, 15.0, horizon
    )
    assert sentinel[0] == explicit[0] == 0
    assert sentinel == explicit
