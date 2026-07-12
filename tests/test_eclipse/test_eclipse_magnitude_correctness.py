"""Correctness of solar-eclipse magnitude and central-path width.

The existing tests range-check these (in places with a 3:1 tolerance). This pins
the magnitude to an independent geometric reimplementation from the library's own
topocentric apparent Sun/Moon (centres' angular separation vs their angular
radii), classifies total vs annular by magnitude, and sanity-bounds the path
width on the central line.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import SUN, MOON, FLG_SWIEPH, FLG_EQUATORIAL, FLG_TOPOCTR

_RSUN_KM = 696000.0
_RMOON_KM = 1738.15
_AU_KM = 149597870.7

# (label, total?, search-from) — two total + one annular.
_ECLIPSES = [
    ("2017-08-21", True, (2017, 8, 1)),
    ("2024-04-08", True, (2024, 4, 1)),
    ("2023-10-14", False, (2023, 10, 1)),  # annular
]


def _ang_sep(ra1, d1, ra2, d2):
    ra1, d1, ra2, d2 = map(math.radians, (ra1, d1, ra2, d2))
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(ra1 - ra2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _greatest(year, month, day):
    jd0 = L.julday(year, month, day, 0.0, L.GREG_CAL)
    _, tret = L.sol_eclipse_when_glob(jd0, FLG_SWIEPH, 0, False)
    jdmax = tret[0]
    _, geopos, _ = L.sol_eclipse_where(jdmax, FLG_SWIEPH)
    return jdmax, geopos[0], geopos[1]


@pytest.mark.parametrize("label,total,start", _ECLIPSES)
def test_eclipse_magnitude_matches_geometry(label, total, start):
    jdmax, lon, lat = _greatest(*start)
    _, attr = L.sol_eclipse_how(jdmax, (lon, lat, 0.0), FLG_SWIEPH)
    lib_mag = attr[0]

    # Independent magnitude from the library's own topocentric Sun/Moon.
    L.set_topo(lon, lat, 0.0)
    base = FLG_SWIEPH | FLG_EQUATORIAL | FLG_TOPOCTR
    s, _ = L.calc_ut(jdmax, SUN, base)
    m, _ = L.calc_ut(jdmax, MOON, base)
    L.set_topo(0.0, 0.0, 0.0)
    r_sun = math.degrees(math.asin(_RSUN_KM / (s[2] * _AU_KM)))
    r_moon = math.degrees(math.asin(_RMOON_KM / (m[2] * _AU_KM)))
    sep = _ang_sep(s[0], s[1], m[0], m[1])
    ind_mag = (r_sun + r_moon - sep) / (2.0 * r_sun)

    assert abs(lib_mag - ind_mag) < 5e-3, f"{label}: {lib_mag} vs {ind_mag}"
    # Total -> magnitude > 1 (Moon covers the whole disc); annular -> < 1.
    assert (lib_mag > 1.0) == total, f"{label}: mag={lib_mag} total={total}"


@pytest.mark.parametrize("label,total,start", _ECLIPSES)
def test_eclipse_path_width_sane(label, total, start):
    jdmax, lon, lat = _greatest(*start)
    width = L.calc_eclipse_path_width(jdmax, lat, lon, FLG_SWIEPH)
    # Central path of a total/annular eclipse: tens to a few hundred km, positive.
    assert 30.0 < width < 500.0, f"{label}: path width {width} km"
