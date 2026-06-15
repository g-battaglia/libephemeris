"""Correctness of the Besselian elements (calc_besselian_*).

The per-element tests (test_besselian_x/y/d/...) only check the return type and a
wide plausible range. This validates the actual numbers two independent ways:

  * an in-test reimplementation of the fundamental-plane geometry from the
    library's own apparent Sun/Moon equatorial positions (a separate code path
    from eclipse.py::_besselian_core, so it catches projection / sign / constant
    bugs there) — must agree to round-off;
  * physical invariants at greatest eclipse (penumbra radius ~0.53-0.56 Earth
    radii; umbra radius l2 < 0 for a total eclipse; axis declination within a few
    degrees of the Moon's declination).
"""
from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import SUN, MOON, FLG_SWIEPH, FLG_EQUATORIAL
from libephemeris.time_utils import sidtime
import libephemeris.eclipse as ecl

_RE = ecl._ECL_REARTH_AU
_RSUN = ecl._ECL_RSUN_AU
_K_PEN = 0.2725076
_K_UMB = 0.272281

# (label, total, search-from) for a few well-separated eclipses.
_ECLIPSES = [
    ("2017-08-21", True, (2017, 8, 1)),
    ("2024-04-08", True, (2024, 4, 1)),
    ("2019-07-02", True, (2019, 7, 1)),
]


def _greatest_eclipse_jd(year, month, day):
    jd0 = L.julday(year, month, day, 0.0, L.GREG_CAL)
    _, tret = L.sol_eclipse_when_glob(jd0, FLG_SWIEPH, 0, False)
    return tret[0]


def _independent_besselian(jd_ut):
    """Reimplement (x, y, d, mu, l1, l2) from the library's apparent equatorial
    Sun/Moon — independent of eclipse.py::_besselian_core."""
    base = FLG_SWIEPH | FLG_EQUATORIAL
    sp, _ = L.calc_ut(jd_ut, SUN, base)
    mp, _ = L.calc_ut(jd_ut, MOON, base)

    def cart(p):
        ra, dec, r = math.radians(p[0]), math.radians(p[1]), p[2]
        return (r * math.cos(dec) * math.cos(ra),
                r * math.cos(dec) * math.sin(ra),
                r * math.sin(dec))

    rs, rm = cart(sp), cart(mp)
    g = [rs[i] - rm[i] for i in range(3)]
    glen = math.sqrt(sum(c * c for c in g))
    gh = [c / glen for c in g]
    d = math.asin(max(-1.0, min(1.0, gh[2])))
    a = math.atan2(gh[1], gh[0])
    exy = math.hypot(gh[0], gh[1])
    ex = [-gh[1] / exy, gh[0] / exy, 0.0]
    ey = [gh[1] * ex[2] - gh[2] * ex[1],
          gh[2] * ex[0] - gh[0] * ex[2],
          gh[0] * ex[1] - gh[1] * ex[0]]
    x = sum(rm[i] * ex[i] for i in range(3)) / _RE
    y = sum(rm[i] * ey[i] for i in range(3)) / _RE
    z = sum(rm[i] * gh[i] for i in range(3)) / _RE
    mu = (sidtime(jd_ut) * 15.0 - math.degrees(a)) % 360.0
    g_er, rsun_er = glen / _RE, _RSUN / _RE
    sf1, sf2 = (rsun_er + _K_PEN) / g_er, (rsun_er - _K_UMB) / g_er
    l1 = z * sf1 / math.sqrt(1 - sf1 * sf1) + _K_PEN / math.sqrt(1 - sf1 * sf1)
    l2 = z * sf2 / math.sqrt(1 - sf2 * sf2) - _K_UMB / math.sqrt(1 - sf2 * sf2)
    return x, y, math.degrees(d), mu, l1, l2


@pytest.mark.parametrize("label,total,start", _ECLIPSES)
def test_besselian_matches_independent_geometry(label, total, start):
    jd = _greatest_eclipse_jd(*start)
    lib = (L.calc_besselian_x(jd), L.calc_besselian_y(jd), L.calc_besselian_d(jd),
           L.calc_besselian_mu(jd), L.calc_besselian_l1(jd), L.calc_besselian_l2(jd))
    ind = _independent_besselian(jd)
    # Same positions, independent geometry -> agree to round-off.
    assert abs(lib[0] - ind[0]) < 1e-6, f"{label} x"
    assert abs(lib[1] - ind[1]) < 1e-6, f"{label} y"
    assert abs(lib[2] - ind[2]) < 1e-6, f"{label} d"
    assert abs((lib[3] - ind[3] + 180) % 360 - 180) < 1e-6, f"{label} mu"
    assert abs(lib[4] - ind[4]) < 1e-6, f"{label} l1"
    assert abs(lib[5] - ind[5]) < 1e-6, f"{label} l2"


@pytest.mark.parametrize("label,total,start", _ECLIPSES)
def test_besselian_physical_invariants(label, total, start):
    jd = _greatest_eclipse_jd(*start)
    l1 = L.calc_besselian_l1(jd)
    l2 = L.calc_besselian_l2(jd)
    d = L.calc_besselian_d(jd)
    moon_dec = L.calc_ut(jd, MOON, FLG_SWIEPH | FLG_EQUATORIAL)[0][1]
    # Penumbral radius on the fundamental plane is ~0.53-0.56 Earth radii.
    assert 0.52 < l1 < 0.57, f"{label}: l1={l1}"
    # A total eclipse has the umbral vertex beyond the plane -> l2 < 0.
    if total:
        assert l2 < 0, f"{label}: total eclipse should have l2<0, got {l2}"
    # The shadow-axis declination tracks the Moon's at greatest eclipse.
    assert abs(d - moon_dec) < 3.0, f"{label}: d={d} vs Moon dec={moon_dec}"
