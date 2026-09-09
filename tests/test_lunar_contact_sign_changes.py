# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Check the physical sign change on either side of lunar contacts."""

from __future__ import annotations

import pytest

from libephemeris import eclipse
from libephemeris.constants import ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL


@pytest.mark.parametrize(
    "jd_max,expected_class",
    [
        (2459715.6746893055, ECL_TOTAL),
        (2459537.877033952, ECL_PARTIAL),
        (2460070.2243118742, ECL_PENUMBRAL),
    ],
)
def test_contact_residuals_change_sign(jd_max, expected_class):
    """Ingress and egress cross the oriented boundary of the relevant cone."""
    word, _attr, _geometry = eclipse._lun_how_core(jd_max)
    assert word == expected_class
    contacts = eclipse._lun_eclipse_phase_times(jd_max, word)
    pairs = [(6, 7, "penumbral")]
    if word & (ECL_TOTAL | ECL_PARTIAL):
        pairs.append((2, 3, "partial"))
    else:
        assert contacts[2] == contacts[3] == 0.0
    if word & ECL_TOTAL:
        pairs.append((4, 5, "total"))
    else:
        assert contacts[4] == contacts[5] == 0.0

    def residual(jd, phase):
        _word, _values, geometry = eclipse._lun_how_core(jd)
        if phase == "penumbral":
            half_width = geometry.penumbral_plane_diameter_km / 2.0
            limb = geometry.shadowed_radius_km / geometry.cos_penumbral_half_angle
        else:
            half_width = -geometry.umbral_plane_diameter_km / 2.0
            limb = geometry.shadowed_radius_km / geometry.cos_umbral_half_angle
            if phase == "total":
                limb = -limb
        return half_width + limb - geometry.axis_offset_km

    # This is temporal sampling, not an output tolerance: both probes are
    # 0.1 seconds from the contact and require strict inequalities.
    step = 0.1 / 86400.0
    for ingress, egress, phase in pairs:
        assert contacts[ingress] < jd_max < contacts[egress]
        assert residual(contacts[ingress] - step, phase) < 0.0
        assert residual(contacts[ingress] + step, phase) > 0.0
        assert residual(contacts[egress] - step, phase) > 0.0
        assert residual(contacts[egress] + step, phase) < 0.0
