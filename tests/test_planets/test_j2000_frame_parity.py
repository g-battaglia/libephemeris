"""Regression tests for the reference J2000 frame convention.

Measured black-box against the reference API (star fit exact to 0.0003"):

* its equatorial J2000 frame is the **mean equator/equinox of J2000**, i.e.
  ICRS rotated by the frame bias (not raw ICRS);
* its ecliptic J2000 frame is that equatorial frame rotated by the **IAU 2006**
  J2000 obliquity (84381.406"), not the IAU 1976 / SPICE-ECLIPJ2000 value
  (84381.448").

An earlier revision used raw ICRS + 84381.448", leaving a constant
~0.036-0.043" latitude offset on every FLG_J2000 ecliptic output (planets,
stars, XYZ) and a matching declination offset on nodes/apsides EQ+J2000; the
osculating inclination was off by exactly 0.042". All values below are frozen
oracle outputs.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe

JD = 2451545.0
ARCSEC = 1.0 / 3600.0


def _dlon(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


@pytest.mark.unit
class TestJ2000FramePlanets:
    """Planet FLG_J2000 outputs against frozen oracle values."""

    def test_mars_ecliptic_j2000(self):
        r, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000)
        # Oracle: 327.967172365, -1.067785441 (fit residual measured 0.0003")
        assert _dlon(r[0], 327.967172365) < 0.01 * ARCSEC
        assert abs(r[1] - (-1.067785441)) < 0.01 * ARCSEC

    def test_mars_equatorial_j2000(self):
        r, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000 | swe.FLG_EQUATORIAL)
        # Oracle: 330.520865902, -13.181931186 (mean equator/equinox J2000,
        # NOT raw ICRS: the frame bias moves Dec here by ~0.02")
        assert _dlon(r[0], 330.520865902) < 0.01 * ARCSEC
        assert abs(r[1] - (-13.181931186)) < 0.01 * ARCSEC

    def test_mars_xyz_j2000(self):
        r, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000 | swe.FLG_XYZ)
        oracle = (1.567790097, -0.980913433, -0.034469475)
        dist = math.sqrt(sum(v * v for v in oracle))
        # Component tolerance equivalent to 0.01" of arc at Mars distance
        tol = dist * 0.01 / 206265.0
        for got, want in zip(r[:3], oracle):
            assert abs(got - want) < tol

    def test_true_node_equatorial_j2000(self):
        r, _ = swe.calc_ut(JD, swe.TRUE_NODE, swe.FLG_J2000 | swe.FLG_EQUATORIAL)
        # Oracle: 126.278989074, 19.264856594. The old raw-ICRS + 84381.448"
        # convention was ~0.034" off in declination here.
        assert _dlon(r[0], 126.278989074) < 0.02 * ARCSEC
        assert abs(r[1] - 19.264856594) < 0.02 * ARCSEC

    def test_ecliptic_equatorial_j2000_consistency(self):
        """Eq J2000 must equal ecl J2000 rotated by exactly -84381.406"."""
        ecl, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000)
        eq, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000 | swe.FLG_EQUATORIAL)
        eps = 84381.406 / 3600.0
        lon_r, lat_r = math.radians(ecl[0]), math.radians(ecl[1])
        e = math.radians(eps)
        x = math.cos(lat_r) * math.cos(lon_r)
        y = math.cos(lat_r) * math.sin(lon_r)
        z = math.sin(lat_r)
        # ecliptic -> equatorial (inverse rotation about x)
        ye = y * math.cos(e) - z * math.sin(e)
        ze = y * math.sin(e) + z * math.cos(e)
        ra = math.degrees(math.atan2(ye, x)) % 360.0
        dec = math.degrees(math.asin(ze))
        assert _dlon(ra, eq[0]) < 0.001 * ARCSEC
        assert abs(dec - eq[1]) < 0.001 * ARCSEC


@pytest.mark.unit
class TestJ2000FrameStars:
    """Fixed-star FLG_J2000 outputs against frozen oracle values."""

    def test_regulus_ecliptic_j2000(self):
        r, _, _ = swe.fixstar2_ut("Regulus", JD, swe.FLG_SWIEPH | swe.FLG_J2000)
        # Oracle: 149.832910673, 0.464813218. The old frame was -0.036" off
        # in latitude (constant across epochs).
        assert _dlon(r[0], 149.832910673) < 0.005 * ARCSEC
        assert abs(r[1] - 0.464813218) < 0.005 * ARCSEC

    def test_vega_ecliptic_j2000(self):
        r, _, _ = swe.fixstar2_ut("Vega", JD, swe.FLG_SWIEPH | swe.FLG_J2000)
        # Oracle: 285.304212269, 61.733258129 (old frame: +0.040" in latitude)
        assert _dlon(r[0], 285.304212269) < 0.01 * ARCSEC
        assert abs(r[1] - 61.733258129) < 0.005 * ARCSEC

    def test_regulus_equatorial_j2000(self):
        r, _, _ = swe.fixstar2_ut(
            "Regulus", JD, swe.FLG_SWIEPH | swe.FLG_J2000 | swe.FLG_EQUATORIAL
        )
        # Oracle: 152.096565050, 11.965851752 (mean equator J2000 = frame
        # bias applied; raw ICRS was ~0.016" off)
        assert _dlon(r[0], 152.096565050) < 0.005 * ARCSEC
        assert abs(r[1] - 11.965851752) < 0.005 * ARCSEC

    def test_j2000_latitude_matches_of_date_at_j2000(self):
        """At the J2000 epoch the J2000-frame latitude equals the of-date
        latitude (the ecliptic planes coincide there); the old frame broke
        this identity by ~0.04"."""
        of_date, _, _ = swe.fixstar2_ut("Regulus", JD, swe.FLG_SWIEPH)
        j2000, _, _ = swe.fixstar2_ut("Regulus", JD, swe.FLG_SWIEPH | swe.FLG_J2000)
        assert abs(of_date[1] - j2000[1]) < 0.005 * ARCSEC


@pytest.mark.unit
class TestJ2000FrameElements:
    """Osculating elements use the reference J2000 ecliptic frame."""

    def test_mars_inclination(self):
        el = swe.get_orbital_elements(JD, swe.MARS, swe.FLG_SWIEPH)
        # Oracle i = 1.8498883; the old 84381.448" frame was off by exactly
        # 0.042" (1.8498765).
        assert abs(el[2] - 1.8498883) < 0.01 * ARCSEC

    def test_mars_ascending_node(self):
        el = swe.get_orbital_elements(JD, swe.MARS, swe.FLG_SWIEPH)
        # Oracle Om = 49.5618579; the old frame was +0.53" off (the 0.042"
        # plane tilt amplified by 1/sin(i)). Residual after the fix: +0.135".
        assert _dlon(el[3], 49.5618579) < 0.3 * ARCSEC
