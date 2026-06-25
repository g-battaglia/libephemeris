"""LEB -> Skyfield fallback *correctness*.

When the LEB fast path can't serve a (body, flag) combination it raises and the
dispatcher silently falls through to Skyfield. The existing tests only assert
that the LEB path raises KeyError — never that the fallback returns the *right*
number. This checks that, for every combination that falls back, the result is
identical to forcing the Skyfield backend directly (the fallback must BE
Skyfield, not an approximation of it).
"""
from __future__ import annotations

import os

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN, MOON, MARS, MEAN_NODE, MEAN_APOG, FLG_SWIEPH, FLG_TOPOCTR, FLG_ICRS,
    FLG_RADIANS,
)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
LEB_MEDIUM = os.path.join(PROJECT_ROOT, "data", "leb", "ephemeris_medium.leb")
CUPIDO = 40
JD = 2451545.0

skip_no_leb = pytest.mark.skipif(
    not os.path.exists(LEB_MEDIUM), reason="medium LEB not found"
)

# (label, body, flags) — each is a combination the LEB path does not serve and
# therefore must route to Skyfield (topocentric ecliptic-direct / heliocentric
# pipelines, ICRS, radians-on-fallback bodies).
_FALLBACK_CASES = [
    ("node_topo", MEAN_NODE, FLG_SWIEPH | FLG_TOPOCTR),
    ("lilith_topo", MEAN_APOG, FLG_SWIEPH | FLG_TOPOCTR),
    ("uranian_topo", CUPIDO, FLG_SWIEPH | FLG_TOPOCTR),
    ("mars_icrs", MARS, FLG_SWIEPH | FLG_ICRS),
    ("moon_icrs", MOON, FLG_SWIEPH | FLG_ICRS),
    ("sun_radians", SUN, FLG_SWIEPH | FLG_RADIANS),
]


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
    assert worst < 1e-9, f"{label}: LEB fallback != Skyfield by {worst:.2e}"
