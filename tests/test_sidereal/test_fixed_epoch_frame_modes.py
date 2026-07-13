"""Independent frame invariants for the standard fixed-epoch modes."""

from __future__ import annotations

from collections.abc import Iterator
import math

import erfa
import numpy as np
import pytest

import libephemeris as le
from libephemeris.precession_vondrak import vondrak_precession_matrix
from libephemeris.sidereal_epoch import FIXED_EPOCH_T0


J2000 = 2451545.0


@pytest.fixture(autouse=True)
def _reset_state() -> Iterator[None]:
    le.reset_session()
    yield
    le.reset_session()


def _rot_x(angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array(((1.0, 0.0, 0.0), (0.0, c, s), (0.0, -s, c)))


def _standard_epoch_matrix(mode: int, equatorial: bool) -> np.ndarray:
    t0 = FIXED_EPOCH_T0[mode]
    p_t0 = np.asarray(vondrak_precession_matrix(t0))
    p_j2000 = np.asarray(vondrak_precession_matrix(J2000))
    equatorial_matrix = p_t0 @ p_j2000.T
    if equatorial:
        return equatorial_matrix
    return (
        _rot_x(float(erfa.obl06(t0, 0.0)))
        @ equatorial_matrix
        @ _rot_x(-float(erfa.obl06(J2000, 0.0)))
    )


def _erfa_epoch_jd(epoch: float, *, besselian: bool = False) -> float:
    parts = erfa.epb2jd(epoch) if besselian else erfa.epj2jd(epoch)
    return float(parts[0] + parts[1])


def test_fixed_reference_epochs_come_from_erfa_definitions() -> None:
    assert FIXED_EPOCH_T0[le.SIDM_J2000] == _erfa_epoch_jd(2000.0)
    assert FIXED_EPOCH_T0[le.SIDM_J1900] == _erfa_epoch_jd(1900.0)
    assert FIXED_EPOCH_T0[le.SIDM_B1950] == _erfa_epoch_jd(1950.0, besselian=True)


@pytest.mark.parametrize("body", [le.MOON, le.JUPITER])
@pytest.mark.parametrize(
    "representation",
    [0, le.FLG_EQUATORIAL, le.FLG_XYZ, le.FLG_RADIANS],
)
def test_j2000_sidereal_mode_is_j2000_mean_frame(
    body: int, representation: int
) -> None:
    common = le.FLG_SPEED | representation
    le.set_sid_mode(le.SIDM_J2000)
    sidereal, sidereal_flags = le.calc_ut(2448000.5, body, common | le.FLG_SIDEREAL)
    tropical, _ = le.calc_ut(
        2448000.5,
        body,
        common | le.FLG_J2000 | le.FLG_NONUT,
    )
    assert sidereal == pytest.approx(tropical, abs=2e-12)
    assert sidereal_flags & le.FLG_SIDEREAL
    assert sidereal_flags & le.FLG_NONUT
    assert not (sidereal_flags & le.FLG_J2000)


@pytest.mark.parametrize("mode", [le.SIDM_J1900, le.SIDM_B1950])
@pytest.mark.parametrize("equatorial", [False, True])
@pytest.mark.parametrize("body", [le.MOON, le.JUPITER])
def test_standard_epoch_xyz_is_published_frame_rotation(
    mode: int, equatorial: bool, body: int
) -> None:
    channel = le.FLG_EQUATORIAL if equatorial else 0
    base, _ = le.calc_ut(
        2448000.5,
        body,
        le.FLG_SPEED | le.FLG_XYZ | le.FLG_J2000 | le.FLG_NONUT | channel,
    )
    le.set_sid_mode(mode)
    actual, _ = le.calc_ut(
        2448000.5,
        body,
        le.FLG_SPEED | le.FLG_XYZ | le.FLG_SIDEREAL | channel,
    )
    matrix = _standard_epoch_matrix(mode, equatorial)
    expected_position = matrix @ np.asarray(base[:3])
    expected_velocity = matrix @ np.asarray(base[3:])
    assert actual[:3] == pytest.approx(expected_position, abs=2e-12)
    assert actual[3:] == pytest.approx(expected_velocity, abs=2e-12)


@pytest.mark.parametrize("mode", [le.SIDM_J1900, le.SIDM_B1950])
@pytest.mark.parametrize("equatorial", [False, True])
def test_spherical_and_xyz_fixed_epoch_channels_are_identical(
    mode: int, equatorial: bool
) -> None:
    channel = le.FLG_EQUATORIAL if equatorial else 0
    le.set_sid_mode(mode)
    spherical, _ = le.calc_ut(
        2460310.5,
        le.SATURN,
        le.FLG_SPEED | le.FLG_SIDEREAL | channel,
    )
    cartesian, _ = le.calc_ut(
        2460310.5,
        le.SATURN,
        le.FLG_SPEED | le.FLG_SIDEREAL | le.FLG_XYZ | channel,
    )
    lon, lat, distance, dlon, dlat, radial_speed = spherical
    lon = math.radians(lon)
    lat = math.radians(lat)
    dlon = math.radians(dlon)
    dlat = math.radians(dlat)
    cl, sl = math.cos(lon), math.sin(lon)
    cb, sb = math.cos(lat), math.sin(lat)
    expected = (
        distance * cb * cl,
        distance * cb * sl,
        distance * sb,
        radial_speed * cb * cl - distance * sb * dlat * cl - distance * cb * sl * dlon,
        radial_speed * cb * sl - distance * sb * dlat * sl + distance * cb * cl * dlon,
        radial_speed * sb + distance * cb * dlat,
    )
    # The spherical channel carries numerically differentiated angular rates;
    # converting them back to Cartesian introduces a few 1e-12 AU/day of
    # round-off relative to the directly rotated Cartesian state.
    assert cartesian == pytest.approx(expected, abs=1e-11)
