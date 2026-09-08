# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Verifica i residui lunari effettivi con geometrie sintetiche esatte."""

from __future__ import annotations

import pytest

from libephemeris import eclipse
from libephemeris.constants import ECL_TOTAL
from libephemeris.shadow_geometry import ShadowGeometry


@pytest.mark.parametrize(
    "offset,inner_diameter,expected",
    [
        (64.0, -32.0, (32.0, -32.0, -64.0)),
        (0.125, -32.0, (95.875, 31.875, -0.125)),
        (8.0, 16.0, (88.0, 0.0, -32.0)),
    ],
)
def test_phase_search_uses_signed_record_fields(
    monkeypatch, offset, inner_diameter, expected
):
    """Esercita i residui reali con ombra, asse vicino e antombra sintetici."""
    geometry = ShadowGeometry(
        offset,
        umbral_plane_diameter_km=inner_diameter,
        penumbral_plane_diameter_km=128.0,
        cos_umbral_half_angle=0.5,
        cos_penumbral_half_angle=0.25,
        shadowed_radius_km=8.0,
        shadow_misses_body=False,
    )
    geometry.validate()
    evaluations = []
    calls = []
    jd_max = 2451545.0
    flags = 2

    def core(jd, actual_flags):
        calls.append((jd, actual_flags))
        return ECL_TOTAL, [], geometry

    def root(residual, lo, hi):
        evaluations.append(residual(jd_max))
        return 0.0

    monkeypatch.setattr(eclipse, "_lun_how_core", core)
    monkeypatch.setattr(eclipse, "_root_bisect", root)
    result = eclipse._lun_eclipse_phase_times(jd_max, ECL_TOTAL, flags)

    assert evaluations == [value for value in expected for _ in range(2)]
    assert calls == [(jd_max, flags)] * 6
    assert result == (jd_max,) + (0.0,) * 9
