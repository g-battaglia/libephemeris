# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""SIDM_GALALIGN_MARDYKS (34) is a fixed-epoch *frame* mode (measured
black-box, see sidereal_epoch.py): every sidereal output is the apparent
position projected onto the mean frames of t0 = JD 2451079.771 (the
September-1998 galactic-alignment equinox), with a constant 30-degree
longitude offset on the ecliptic/XYZ channels and no offset on the
equatorial channel.  A scalar ayanamsha subtraction reproduces the
longitudes only near J2000 and leaves the latitudes on the wrong plane
(~24" off at 1946 from the precession of the ecliptic pole).

Reference values captured from the reference API (mode 34, various epochs).
"""

from __future__ import annotations

import pytest

import libephemeris as le

SIDEQ = le.FLG_SIDEREAL | le.FLG_EQUATORIAL | le.FLG_SPEED
SIDECL = le.FLG_SIDEREAL | le.FLG_SPEED


@pytest.fixture(autouse=True)
def _galalign_mode():
    le.set_sid_mode(le.SIDM_GALALIGN_MARDYKS, 0, 0)
    yield
    le.set_sid_mode(0)


def _assert_lonlat(got, want_lon, want_lat, tol_arcsec):
    dlon = abs((got[0] - want_lon + 180.0) % 360.0 - 180.0) * 3600.0
    dlat = abs(got[1] - want_lat) * 3600.0
    assert dlon < tol_arcsec, (got[0], want_lon, dlon)
    assert dlat < tol_arcsec, (got[1], want_lat, dlat)


# (body, jd_ut, channel) -> (lon/RA, lat/Dec) from the reference API.
ECLIPTIC_CASES = [
    ("SUN", 2432000.5, 67.38287177082623, 0.006650248019015653, 0.05),
    ("MARS", 2432000.5, 125.68183348515672, 1.0632764787030582, 0.05),
    ("MEAN_NODE", 2432000.5, 50.72960896486768, 0.006803320970692691, 0.1),
    ("TRUE_NODE", 2432000.5, 51.526986748618796, 0.006809712470273666, 0.05),
    ("MEAN_APOG", 2432000.5, 216.80798570248004, 1.234620405386485, 0.25),
    ("SATURN", 2451545.0, 10.381734482922692, -2.4449733901716972, 0.1),
    ("SUN", 2466000.5, 97.12578428296837, -0.003796817438453902, 0.15),
]

EQUATORIAL_CASES = [
    ("MARS", 2451545.0, 330.5037285, -13.1881043, 0.05),
    ("MEAN_NODE", 2451545.0, 127.3777185, 19.0100205, 0.05),
    ("MEAN_APOG", 2451545.0, 263.0466469, -19.8616131, 0.4),
    ("TRUE_NODE", 2451545.0, 126.2606709, 19.269052, 0.05),
]


@pytest.mark.parametrize("name,jd,lon,lat,tol", ECLIPTIC_CASES)
def test_ecliptic_channel_projects_onto_t0_plane(name, jd, lon, lat, tol):
    got = le.calc_ut(jd, getattr(le, name), SIDECL)[0]
    _assert_lonlat(got, lon, lat, tol)


@pytest.mark.parametrize("name,jd,ra,dec,tol", EQUATORIAL_CASES)
def test_equatorial_channel_is_t0_mean_equator(name, jd, ra, dec, tol):
    got = le.calc_ut(jd, getattr(le, name), SIDEQ)[0]
    _assert_lonlat(got, ra, dec, tol)


def test_sun_latitude_nonzero_off_t0_plane():
    """The tell-tale of the frame model: the Sun's sidereal latitude is the
    tilt between the ecliptic of date and the ecliptic of t0 (~24" in 1946),
    not zero as a scalar-ayanamsha model would give."""
    lat = le.calc_ut(2432000.5, le.SUN, SIDECL)[0][1]
    assert abs(lat - 0.0066502) * 3600.0 < 0.05
    assert lat > 0.005  # far from the scalar model's ~0


def test_xyz_channels_consistent_with_spherical():
    import math

    for flags in (SIDECL, SIDEQ):
        sph = le.calc_ut(2432000.5, le.SATURN, flags)[0]
        xyz = le.calc_ut(2432000.5, le.SATURN, flags | le.FLG_XYZ)[0]
        lon = math.degrees(math.atan2(xyz[1], xyz[0])) % 360.0
        r = math.sqrt(sum(c * c for c in xyz[:3]))
        lat = math.degrees(math.asin(xyz[2] / r))
        assert abs((lon - sph[0] + 180.0) % 360.0 - 180.0) < 1e-9
        assert abs(lat - sph[1]) < 1e-9
        assert abs(r - sph[2]) < 1e-12


def test_retflag_echoes_nonut_like_fixed_epoch_modes():
    _, rf = le.calc_ut(2451545.0, le.SUN, SIDECL)
    assert rf & le.FLG_NONUT
    assert rf & le.FLG_SIDEREAL
    assert not (rf & le.FLG_J2000)


def test_scalar_ayanamsha_value_unchanged():
    """get_ayanamsa keeps returning the mode's scalar value (30 deg at t0
    plus accumulated precession); the frame model does not remove it."""
    assert abs(le.get_ayanamsa_ut(2451545.0) - 30.01779390499155) < 1e-4


STAR_CASES = [
    ("Regulus", 2432000.5, SIDECL, 119.8114143, 0.4660363),
    ("Regulus", 2432000.5, SIDEQ, 152.0764644, 11.9746303),
    ("Spica", 2466000.5, SIDECL, 173.8219281, -2.055088),
    ("Spica", 2466000.5, SIDEQ, 201.2796126, -11.1547307),
]


@pytest.mark.parametrize("star,jd,flags,lon,lat", STAR_CASES)
def test_fixed_stars_follow_the_frame_model(star, jd, flags, lon, lat):
    got = le.fixstar2_ut(star, jd, flags)[0]
    _assert_lonlat(got, lon, lat, 0.05)


HOUSE_CASES = [
    # (jd, lat, lon, hsys, cusp1, mc)
    (2432000.5, 48.0, 11.0, "P", 5.446569, 256.841836),
    (2451545.0, 48.0, 11.0, "P", 11.576649, 259.816188),
    (2466000.5, -35.0, 151.0, "K", 164.168975, 68.145885),
    (2451545.0, 65.0, -20.0, "I", 257.340526, 231.217744),
    # Near-t0 epochs bracket the sign-quirk window of the fixed-epoch path.
    (2451079.0, 48.0, 11.0, "P", 225.7301649310638, 163.20470506989548),
    (2451081.0, 48.0, 11.0, "P", 227.33140084210652, 165.32743033879328),
]


@pytest.mark.parametrize("jd,glat,glon,hsys,cusp1,mc", HOUSE_CASES)
def test_houses_follow_the_frame_model(jd, glat, glon, hsys, cusp1, mc):
    cusps, ascmc = le.houses_ex(jd, glat, glon, ord(hsys), le.FLG_SIDEREAL)
    assert abs((cusps[0] - cusp1 + 180.0) % 360.0 - 180.0) * 3600.0 < 0.2
    assert abs((ascmc[1] - mc + 180.0) % 360.0 - 180.0) * 3600.0 < 0.2


def test_sunshine_makransky_stays_distinct_in_fixed_epoch_mode():
    """Fixed-epoch modes keep 'i' (Makransky) as a solution distinct from
    'I' (Treindl) — the ayanamsha-based sidereal path collapses them, the
    fixed-epoch path must not (reference values, cusps 2/6 at lat 65)."""
    cusps, _ = le.houses_ex(2451545.0, 65.0, -20.0, ord("i"), le.FLG_SIDEREAL)
    assert abs(cusps[1] - 294.2779641928922) * 3600.0 < 0.5
    assert abs(cusps[5] - 163.12088902642247) * 3600.0 < 0.5
    c_i, _ = le.houses_ex(2451545.0, 65.0, -20.0, ord("I"), le.FLG_SIDEREAL)
    assert abs((cusps[1] - c_i[1] + 180.0) % 360.0 - 180.0) > 90.0


def test_leb_backend_matches_skyfield_backend():
    for flags in (SIDECL, SIDEQ, SIDECL | le.FLG_XYZ):
        le.set_calc_mode("skyfield")
        a = le.calc_ut(2451545.0, le.MARS, flags)[0]
        le.set_calc_mode("leb")
        b = le.calc_ut(2451545.0, le.MARS, flags)[0]
        le.set_calc_mode("auto")
        for i in range(6):
            assert abs(a[i] - b[i]) < 5e-6, (flags, i, a[i], b[i])


def test_context_api_honours_frame_mode():
    from libephemeris import EphemerisContext

    ctx = EphemerisContext()
    ctx.set_sid_mode(le.SIDM_GALALIGN_MARDYKS, 0, 0)
    got = ctx.calc_ut(2432000.5, le.SUN, SIDECL)[0]
    _assert_lonlat(got, 67.38287177082623, 0.006650248019015653, 0.05)
