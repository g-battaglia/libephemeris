"""Fixed-star catalog integrity + fixstar/fixstar2 agreement on named stars.

A full-catalog sweep confirmed every entry returns finite, in-range coordinates,
and that fixstar (legacy) and fixstar2 (modern) agree exactly for proper/Bayer
names. They diverge only for fragile Flamsteed designations (documented in
docs/comparison/known-differences.md §4, Fixed stars); these guards lock the parts that must hold.
"""
from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import FLG_SWIEPH

J2000 = 2451545.0
_NAMED = ["Sirius", "Aldebaran", "Regulus", "Spica", "Antares", "Betelgeuse",
          "Rigel", "Vega", "Polaris", "Arcturus", "Canopus", "Procyon"]


@pytest.fixture(autouse=True)
def _skyfield():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")


@pytest.mark.parametrize("name", _NAMED)
def test_fixstar_v1_v2_agree_for_named_stars(name):
    c1 = L.fixstar_ut(name, J2000, FLG_SWIEPH)[0]
    c2 = L.fixstar2_ut(name, J2000, FLG_SWIEPH)[0]
    d = abs(c1[0] - c2[0]) % 360.0
    assert min(d, 360.0 - d) * 3600.0 < 1.0
    assert abs(c1[1] - c2[1]) * 3600.0 < 1.0


@pytest.mark.slow
def test_full_catalog_no_anomalies():
    """Every fixstar2 entry returns finite, in-range coordinates and dist > 0."""
    stars = L.list_fixed_stars()
    assert len(stars) > 1000
    bad = []
    for s in stars:
        name = s.name if hasattr(s, "name") else str(s)
        comps = L.fixstar2_ut(name, J2000, FLG_SWIEPH)[0]
        if not all(math.isfinite(v) for v in comps):
            bad.append((name, "nan"))
        elif not (0.0 <= comps[0] < 360.001 and -90.001 <= comps[1] <= 90.001
                  and comps[2] > 0.0):
            bad.append((name, comps[:3]))
    assert not bad, f"{len(bad)} catalog anomalies, e.g. {bad[:5]}"
