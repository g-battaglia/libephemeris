"""Independent consistency tests for the IAU J2000 frame.

No expected ephemeris output is frozen here. The tests use the IAU 2006 J2000
obliquity and coordinate identities to cross-check public representations.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe


JD = 2451545.0
IAU_2006_EPS_J2000 = math.radians(84381.406 / 3600.0)


def _angular_difference(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def _ecliptic_to_equatorial(lon: float, lat: float) -> tuple[float, float]:
    lon_r = math.radians(lon)
    lat_r = math.radians(lat)
    x = math.cos(lat_r) * math.cos(lon_r)
    y = math.cos(lat_r) * math.sin(lon_r)
    z = math.sin(lat_r)
    ye = y * math.cos(IAU_2006_EPS_J2000) - z * math.sin(IAU_2006_EPS_J2000)
    ze = y * math.sin(IAU_2006_EPS_J2000) + z * math.cos(IAU_2006_EPS_J2000)
    return math.degrees(math.atan2(ye, x)) % 360.0, math.degrees(
        math.atan2(ze, math.hypot(x, ye))
    )


@pytest.mark.unit
@pytest.mark.parametrize("body", [swe.MARS, swe.TRUE_NODE])
def test_j2000_ecliptic_and_equatorial_are_iau_rotations(body: int) -> None:
    ecl, _ = swe.calc_ut(JD, body, swe.FLG_J2000)
    eq, _ = swe.calc_ut(JD, body, swe.FLG_J2000 | swe.FLG_EQUATORIAL)
    ra, dec = _ecliptic_to_equatorial(ecl[0], ecl[1])
    assert _angular_difference(ra, eq[0]) < 1e-7
    assert abs(dec - eq[1]) < 1e-7


@pytest.mark.unit
def test_j2000_xyz_and_spherical_planet_outputs_are_identical_vectors() -> None:
    spherical, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000)
    xyz, _ = swe.calc_ut(JD, swe.MARS, swe.FLG_J2000 | swe.FLG_XYZ)
    lon = math.radians(spherical[0])
    lat = math.radians(spherical[1])
    radius = spherical[2]
    expected = (
        radius * math.cos(lat) * math.cos(lon),
        radius * math.cos(lat) * math.sin(lon),
        radius * math.sin(lat),
    )
    assert xyz[:3] == pytest.approx(expected, rel=0.0, abs=2e-12)


@pytest.mark.unit
@pytest.mark.parametrize("name", ["Regulus", "Vega"])
def test_j2000_fixed_star_frames_are_self_consistent(name: str) -> None:
    ecl, _, _ = swe.fixstar2_ut(name, JD, swe.FLG_SWIEPH | swe.FLG_J2000)
    eq, _, _ = swe.fixstar2_ut(
        name, JD, swe.FLG_SWIEPH | swe.FLG_J2000 | swe.FLG_EQUATORIAL
    )
    ra, dec = _ecliptic_to_equatorial(ecl[0], ecl[1])
    assert _angular_difference(ra, eq[0]) < 1e-7
    assert abs(dec - eq[1]) < 1e-7


@pytest.mark.unit
def test_osculating_elements_have_physical_ranges() -> None:
    elements = swe.get_orbital_elements(JD, swe.MARS, swe.FLG_SWIEPH)
    assert elements[0] > 0.0
    assert 0.0 <= elements[1] < 1.0
    assert 0.0 <= elements[2] <= 180.0
    assert 0.0 <= elements[3] < 360.0
