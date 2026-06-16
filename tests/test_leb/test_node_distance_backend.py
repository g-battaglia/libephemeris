"""True Node: angular position correct on every backend; distance is Skyfield-only.

A flag-pair differential (LEB vs Skyfield, FLG_XYZ) surfaced that the True Node's
distance — and thus its cartesian XYZ — differs between backends: the Skyfield/
default path computes the osculating node radius (matching the reference ephemeris), while the
LEB path returns a stored proxy (documented limitation, divergences.md 7). The
node's longitude/latitude are correct everywhere. These guards lock that in:
the angular position must agree across backends, and the Skyfield-mode distance
must be the osculating-conic value (clearly different from the LEB proxy).
"""
from __future__ import annotations

import os

import pytest

import libephemeris as L
from libephemeris.constants import TRUE_NODE, MEAN_NODE, OSCU_APOG, FLG_SWIEPH

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
LEB_MEDIUM = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_medium.leb")

skip_no_leb = pytest.mark.skipif(
    not os.path.exists(LEB_MEDIUM), reason="medium LEB not found"
)
_DATES = [2451545.0 + (y - 2000) * 365.25 for y in (1980, 2000, 2020)]


def _wrap(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


@skip_no_leb
@pytest.mark.parametrize("body", [TRUE_NODE, MEAN_NODE, OSCU_APOG])
def test_node_longitude_agrees_across_backends(body):
    """Longitude/latitude must match LEB vs Skyfield for every lunar point."""
    L.set_leb_file(LEB_MEDIUM)
    try:
        for jd in _DATES:
            L.set_calc_mode("leb")
            leb = L.calc_ut(jd, body, FLG_SWIEPH)[0]
            L.set_calc_mode("skyfield")
            sky = L.calc_ut(jd, body, FLG_SWIEPH)[0]
            assert _wrap(leb[0], sky[0]) * 3600.0 < 0.5
            assert abs(leb[1] - sky[1]) * 3600.0 < 0.5
    finally:
        L.set_calc_mode("auto")


def test_true_node_distance_is_osculating_on_skyfield():
    """Default/Skyfield True Node distance is the osculating radius (~0.0024-0.0027
    AU), not the Moon's mean distance — a sanity bound that the correct path holds."""
    L.set_calc_mode("skyfield")
    try:
        for jd in _DATES:
            dist = L.calc_ut(jd, TRUE_NODE, FLG_SWIEPH)[0][2]
            assert 0.0020 < dist < 0.0030, f"node dist {dist}"
    finally:
        L.set_calc_mode("auto")
