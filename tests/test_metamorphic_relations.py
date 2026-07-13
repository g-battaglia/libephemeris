"""Metamorphic relations — self-checking identities that must hold across the API.

These need no external truth: an internal inconsistency between two ways of asking
for the same quantity is itself a bug. They cover the frame/flag conversions
(XYZ, radians, equatorial, helio/geo, south node, sidereal, ET/UT) where a sign
or convention slip would otherwise pass silently.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    MARS,
    JUPITER,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    OSCU_APOG,
    ECL_NUT,
    FLG_SWIEPH,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_XYZ,
    FLG_RADIANS,
    FLG_HELCTR,
    FLG_J2000,
    FLG_TRUEPOS,
)

_BODIES = [SUN, MOON, MERCURY, MARS, JUPITER, PLUTO, MEAN_NODE, TRUE_NODE, OSCU_APOG]
_DATES = [2451545.0 + (y - 2000) * 365.25 for y in (1950, 2000, 2025, 2080)]


def _wrap(a, b):
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


@pytest.fixture(autouse=True)
def _skyfield():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")
    L.set_topo(0.0, 0.0, 0.0)


@pytest.mark.parametrize("body", _BODIES)
def test_xyz_magnitude_equals_distance(body):
    for jd in _DATES:
        sph = L.calc(jd, body, FLG_SWIEPH)[0]
        xyz = L.calc(jd, body, FLG_SWIEPH | FLG_XYZ)[0]
        assert abs(math.hypot(xyz[0], xyz[1], xyz[2]) - sph[2]) < 1e-9


@pytest.mark.parametrize("body", _BODIES)
def test_radians_equals_degrees(body):
    for jd in _DATES:
        d = L.calc(jd, body, FLG_SWIEPH)[0]
        r = L.calc(jd, body, FLG_SWIEPH | FLG_RADIANS)[0]
        assert abs(math.radians(d[0]) - r[0]) < 1e-12
        assert abs(math.radians(d[1]) - r[1]) < 1e-12
        assert abs(d[2] - r[2]) < 1e-12  # distance unchanged


@pytest.mark.parametrize("body", _BODIES)
def test_speed_flag_leaves_position_identical(body):
    for jd in _DATES:
        p0 = L.calc(jd, body, FLG_SWIEPH)[0]
        p1 = L.calc(jd, body, FLG_SWIEPH | FLG_SPEED)[0]
        assert _wrap(p0[0], p1[0]) < 1e-9
        assert abs(p0[1] - p1[1]) < 1e-9 and abs(p0[2] - p1[2]) < 1e-9


@pytest.mark.parametrize("body", [MEAN_NODE, TRUE_NODE])
def test_south_node_is_north_plus_180(body):
    for jd in _DATES:
        n = L.calc(jd, body, FLG_SWIEPH)[0]
        s = L.calc(jd, -body, FLG_SWIEPH)[0]
        assert _wrap((n[0] + 180.0) % 360.0, s[0]) < 1e-6
        assert abs((-n[1]) - s[1]) < 1e-6


@pytest.mark.parametrize("body", [MERCURY, MARS, JUPITER, PLUTO])
def test_helio_plus_sun_equals_geo(body):
    """geocentric = heliocentric + (geocentric Sun), geometrically (FLG_TRUEPOS)."""
    base = FLG_SWIEPH | FLG_J2000 | FLG_XYZ | FLG_TRUEPOS
    for jd in _DATES:
        helio = L.calc(jd, body, base | FLG_HELCTR)[0]
        geo = L.calc(jd, body, base)[0]
        sun = L.calc(jd, SUN, base)[0]
        for i in range(3):
            assert abs((helio[i] + sun[i]) - geo[i]) < 1e-10


@pytest.mark.parametrize("body", _BODIES)
def test_equatorial_equals_rotated_ecliptic(body):
    for jd in _DATES:
        ecl = L.calc(jd, body, FLG_SWIEPH)[0]
        equ = L.calc(jd, body, FLG_SWIEPH | FLG_EQUATORIAL)[0]
        eps = math.radians(L.calc(jd, ECL_NUT, FLG_SWIEPH)[0][0])  # true obliquity
        lam, bet = math.radians(ecl[0]), math.radians(ecl[1])
        ra = math.atan2(
            math.sin(lam) * math.cos(eps) - math.tan(bet) * math.sin(eps), math.cos(lam)
        )
        dec = math.asin(
            math.sin(bet) * math.cos(eps)
            + math.cos(bet) * math.sin(eps) * math.sin(lam)
        )
        assert _wrap(math.degrees(ra) % 360.0, equ[0]) * 3600.0 < 0.05
        assert abs(math.degrees(dec) - equ[1]) * 3600.0 < 0.05


@pytest.mark.parametrize("body", _BODIES)
def test_et_ut_consistency(body):
    """calc(jd_tt) == calc_ut(jd_ut) when jd_tt = jd_ut + Delta-T."""
    for jd in _DATES:
        a = L.calc_ut(jd, body, FLG_SWIEPH)[0]
        c = L.calc(jd + L.deltat(jd), body, FLG_SWIEPH)[0]
        assert _wrap(a[0], c[0]) * 3600.0 < 1e-2
