"""LEB -> Skyfield fallback *correctness*.

When the LEB fast path can't serve a (body, flag) combination it raises and the
dispatcher silently falls through to Skyfield. The existing tests only assert
that the LEB path raises KeyError — never that the fallback returns the *right*
number. This checks that, for every combination that falls back, the result
matches forcing the Skyfield backend directly to within the LEB compression
spec (see the tolerance note below).
"""

from __future__ import annotations

import os

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    MEAN_NODE,
    MEAN_APOG,
    HARRINGTON,
    FLG_SWIEPH,
    FLG_TOPOCTR,
    FLG_ICRS,
    FLG_RADIANS,
)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
LEB_MEDIUM = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_medium.leb")
JD = 2451545.0

skip_no_leb = pytest.mark.skipif(
    not os.path.exists(LEB_MEDIUM), reason="medium LEB not found"
)

# (label, body, flags) — each is a combination the LEB fast path does not serve
# and therefore must route to the Skyfield pipeline (topocentric
# ecliptic-direct / heliocentric pipelines, ICRS, radians-on-fallback bodies).
_FALLBACK_CASES = [
    ("node_topo", MEAN_NODE, FLG_SWIEPH | FLG_TOPOCTR),
    ("lilith_topo", MEAN_APOG, FLG_SWIEPH | FLG_TOPOCTR),
    ("harrington_topo", HARRINGTON, FLG_SWIEPH | FLG_TOPOCTR),
    ("mars_icrs", MARS, FLG_SWIEPH | FLG_ICRS),
    ("moon_icrs", MOON, FLG_SWIEPH | FLG_ICRS),
    ("sun_radians", SUN, FLG_SWIEPH | FLG_RADIANS),
]

# Tolerance: the sealed-mode fallback is NOT bit-identical Skyfield anymore.
# In calculation mode "leb" the dispatcher's fallback runs the Skyfield
# frame-reduction/light-time code, but its persisted vector states come from
# the LEB vector adapter (state._get_computation_ephemeris routes to
# leb_vector.LEBVectorEphemeris) — never from the sealed DE kernel. Every
# fallback case therefore carries the LEB compression error, so the honest
# contract is the LEB spec: < 0.001 arcsec (1 mas) per stored state,
# i.e. 1e-3 / 3600 deg. The Moon shows the largest angular error (~0.16 mas
# at J2000 — same km-level compression error, smallest distance); the other
# cases sit around 1e-13..1e-10 deg, well inside the same spec.
_LEB_SPEC_DEG = 1e-3 / 3600.0  # ~2.78e-7 deg


@skip_no_leb
@pytest.mark.parametrize("label,body,flags", _FALLBACK_CASES)
def test_leb_fallback_equals_skyfield(label, body, flags):
    L.set_leb_file(LEB_MEDIUM)
    L.set_topo(9.0, 45.0, 100.0)
    try:
        L.set_calc_mode("skyfield")
        sky, _ = L.calc_ut(JD, body, flags)
        L.set_calc_mode("leb")  # LEB can't serve -> must fall back to Skyfield
        leb, _ = L.calc_ut(JD, body, flags)
    finally:
        L.set_calc_mode("auto")
        L.set_topo(0.0, 0.0, 0.0)
    worst = max(abs(s - lv) for s, lv in zip(sky[:3], leb[:3]))
    assert worst < _LEB_SPEC_DEG, (
        f"{label}: LEB fallback deviates from Skyfield by {worst:.2e} deg "
        f"(> compression spec {_LEB_SPEC_DEG:.2e} deg)"
    )
