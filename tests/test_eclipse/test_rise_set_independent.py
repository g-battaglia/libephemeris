"""Independent (erfa) validation of rise/set/twilight.

The existing rise/set tests check ordering + loose bounds vs Swiss Ephemeris.
These adjudicate rise_trans against an independent observed-altitude reference:

  * Sun rise + astronomical twilight via erfa.atco13 (observed zenith distance),
    with refraction OFF so lib and erfa share the same altitude condition — pure
    geometry, agree to a few seconds;
  * Moon rise via an independent hour-angle altitude from the library's
    (separately validated) topocentric apparent RA/Dec + erfa sidereal time —
    exercises the ~1 deg diurnal parallax;
  * invariants: rise < transit < set within a day; circumpolar -> res = -2.

Network-free but needs pyerfa; marked slow (root-finding over many positions).
"""

from __future__ import annotations

import math

import pytest

import libephemeris as L
from libephemeris.constants import (
    SUN,
    MOON,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_TOPOCTR,
    CALC_RISE,
    CALC_SET,
    CALC_MTRANSIT,
    BIT_DISC_CENTER,
    BIT_NO_REFRACTION,
    BIT_ASTRO_TWILIGHT,
)

erfa = pytest.importorskip("erfa")

_RMOON_KM = 1738.15
_AU_KM = 149597870.7


@pytest.fixture(autouse=True)
def _skyfield():
    L.set_calc_mode("skyfield")
    yield
    L.set_calc_mode("auto")
    L.set_topo(0.0, 0.0, 0.0)


def _sun_obs_alt(jd_ut, lon, lat, phpa):
    """Sun observed altitude (deg) via erfa.atco13 fed the Sun's ICRS astrometric
    place. phpa=0 disables refraction."""
    jd_tt = jd_ut + L.deltat(jd_ut)
    rd, _ = L.calc(
        jd_tt, SUN, FLG_SWIEPH | FLG_J2000 | FLG_EQUATORIAL | FLG_NOABERR | FLG_NOGDEFL
    )
    rc, dc = math.radians(rd[0]), math.radians(rd[1])
    _aob, zob, _h, _d, _r, _eo = erfa.atco13(
        rc,
        dc,
        0,
        0,
        0,
        0,
        jd_ut,
        0.0,
        0.0,
        math.radians(lon),
        math.radians(lat),
        0.0,
        0,
        0,
        phpa,
        10.0,
        0.5,
        0.55,
    )
    return 90.0 - math.degrees(zob)


def _moon_obs_alt(jd_ut, lon, lat):
    """Moon observed altitude (deg) from its topocentric apparent RA/Dec (lib) and
    erfa sidereal time — independent of lib's rise/set geometry, with parallax."""
    L.set_topo(lon, lat, 0.0)
    rd, _ = L.calc_ut(jd_ut, MOON, FLG_SWIEPH | FLG_EQUATORIAL | FLG_TOPOCTR)
    ra, dec = math.radians(rd[0]), math.radians(rd[1])
    gst = erfa.gst06a(2400000.5, jd_ut - 2400000.5, 2400000.5, jd_ut - 2400000.5)
    H = gst + math.radians(lon) - ra
    phi = math.radians(lat)
    return math.degrees(
        math.asin(
            math.sin(dec) * math.sin(phi) + math.cos(dec) * math.cos(phi) * math.cos(H)
        )
    )


def _find_rise(alt_fn, jd0, target_alt):
    step = 1 / 48.0
    prev = alt_fn(jd0) - target_alt
    jd = jd0
    for _ in range(48 * 2):
        jd += step
        cur = alt_fn(jd) - target_alt
        if prev < 0 <= cur:
            a, b = jd - step, jd
            for _ in range(40):
                m = (a + b) / 2
                if (alt_fn(a) - target_alt < 0) != (alt_fn(m) - target_alt < 0):
                    b = m
                else:
                    a = m
            return (a + b) / 2
        prev = cur
    return None


@pytest.mark.slow
@pytest.mark.parametrize("lon,lat", [(12.5, 41.9), (0.0, 0.0), (-78.5, -0.2)])
def test_sun_rise_and_twilight_match_erfa(lon, lat):
    """Sun rise (disc-center, no refraction) and astronomical twilight match the
    erfa observed-altitude reference to a few seconds (shared geometry)."""
    jd0 = L.julday(2023, 6, 15, 0.0, L.GREG_CAL)
    geopos = (lon, lat, 0.0)

    rf, tret = L.rise_trans(
        jd0, SUN, CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION, geopos
    )
    ref = _find_rise(lambda j: _sun_obs_alt(j, lon, lat, 0.0), jd0, 0.0)
    assert rf == 0 and ref is not None
    assert abs(tret[0] - ref) * 86400.0 < 5.0

    rf2, tret2 = L.rise_trans(
        jd0, SUN, CALC_RISE | BIT_ASTRO_TWILIGHT | BIT_DISC_CENTER, geopos
    )
    ref2 = _find_rise(lambda j: _sun_obs_alt(j, lon, lat, 0.0), jd0, -18.0)
    assert rf2 == 0 and ref2 is not None
    assert abs(tret2[0] - ref2) * 86400.0 < 5.0


@pytest.mark.slow
@pytest.mark.parametrize("lon,lat", [(12.5, 41.9), (0.0, 0.0), (151.2, -33.9)])
def test_moon_rise_matches_independent_hour_angle(lon, lat):
    """Moon rise (disc-center, no refraction) matches the independent hour-angle
    altitude (with ~1 deg diurnal parallax) to a few seconds."""
    jd0 = L.julday(2023, 3, 10, 0.0, L.GREG_CAL)
    geopos = (lon, lat, 0.0)
    rf, tret = L.rise_trans(
        jd0, MOON, CALC_RISE | BIT_DISC_CENTER | BIT_NO_REFRACTION, geopos
    )
    ref = _find_rise(lambda j: _moon_obs_alt(j, lon, lat), jd0, 0.0)
    assert rf == 0 and ref is not None
    assert abs(tret[0] - ref) * 86400.0 < 5.0


def test_rise_transit_set_ordering_and_circumpolar():
    """Invariants: rise < transit < set within a day; a circumpolar Sun returns
    res = -2 (no rise/set) at a polar site in local summer."""
    jd0 = L.julday(2023, 6, 21, 0.0, L.GREG_CAL)
    rr, tr = L.rise_trans(jd0, SUN, CALC_RISE, (12.5, 41.9, 0.0))
    rt, tt = L.rise_trans(jd0, SUN, CALC_MTRANSIT, (12.5, 41.9, 0.0))
    rs, ts = L.rise_trans(jd0, SUN, CALC_SET, (12.5, 41.9, 0.0))
    assert rr == 0 and rt == 0 and rs == 0
    assert tr[0] < tt[0] < ts[0] < jd0 + 1.0
    # Svalbard (78N) in June: midnight sun -> Sun never sets.
    rc, _tc = L.rise_trans(jd0, SUN, CALC_SET, (15.6, 78.2, 0.0))
    assert rc == -2
