"""Independent geometric and representation invariants for ``nod_aps``.

These tests derive expectations from coordinate identities, NASA JPL/Skyfield
observer states, and the Vondrák/IAU frame model used by LibEphemeris.  They do
not contain recorded output from another ephemeris implementation.
"""

from __future__ import annotations

from collections.abc import Iterator
import math

import erfa
import pytest

import libephemeris as eph
from libephemeris.astrometry import apply_aberration_to_position
from libephemeris.precession_vondrak import (
    vondrak_mean_obliquity_rad,
    vondrak_pn_matrix,
    vondrak_precession_matrix,
)
from libephemeris.state import get_planets, get_timescale, get_topo


J2000 = 2451545.0


@pytest.fixture(autouse=True)
def _isolate_observer_and_sidereal_state() -> Iterator[None]:
    eph.reset_session()
    # These invariants dissect the Skyfield reduction internals (get_planets,
    # aberration helpers), so they must not inherit a sealed-LEB mode from the
    # environment (e.g. a local .env exporting LIBEPHEMERIS_MODE=leb).
    eph.set_calc_mode("skyfield")
    yield
    eph.reset_session()


def _spherical_state_to_cartesian(
    state: tuple[float, float, float, float, float, float],
) -> tuple[float, float, float, float, float, float]:
    lon, lat, distance, dlon, dlat, radial_speed = state
    lon = math.radians(lon)
    lat = math.radians(lat)
    dlon = math.radians(dlon)
    dlat = math.radians(dlat)
    cl, sl = math.cos(lon), math.sin(lon)
    cb, sb = math.cos(lat), math.sin(lat)
    return (
        distance * cb * cl,
        distance * cb * sl,
        distance * sb,
        radial_speed * cb * cl - distance * sb * dlat * cl - distance * cb * sl * dlon,
        radial_speed * cb * sl - distance * sb * dlat * sl + distance * cb * cl * dlon,
        radial_speed * sb + distance * cb * dlat,
    )


def _icrs_to_ecliptic(
    vector: tuple[float, float, float], jd_tt: float, *, nutation: bool
) -> tuple[float, float, float]:
    if nutation:
        dpsi, deps = erfa.nut06a(J2000, jd_tt - J2000)
        matrix, obliquity = vondrak_pn_matrix(jd_tt, float(dpsi), float(deps))
    else:
        matrix = vondrak_precession_matrix(jd_tt)
        obliquity = vondrak_mean_obliquity_rad(jd_tt)
    equatorial = tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )
    ce, se = math.cos(obliquity), math.sin(obliquity)
    return (
        equatorial[0],
        equatorial[1] * ce + equatorial[2] * se,
        -equatorial[1] * se + equatorial[2] * ce,
    )


def _ecliptic_date_to_j2000(
    vector: tuple[float, float, float], jd_tt: float
) -> tuple[float, float, float]:
    """Apply the published Vondrák frame matrices as a linear rotation."""
    eps = vondrak_mean_obliquity_rad(J2000)
    ce, se = math.cos(eps), math.sin(eps)
    p_date = vondrak_precession_matrix(jd_tt)
    p_j2000 = vondrak_precession_matrix(J2000)
    equatorial_date = (
        vector[0],
        vector[1] * ce - vector[2] * se,
        vector[1] * se + vector[2] * ce,
    )
    icrs = tuple(
        sum(p_date[row][column] * equatorial_date[row] for row in range(3))
        for column in range(3)
    )
    equatorial_j2000 = tuple(
        sum(p_j2000[row][column] * icrs[column] for column in range(3))
        for row in range(3)
    )
    return (
        equatorial_j2000[0],
        equatorial_j2000[1] * ce + equatorial_j2000[2] * se,
        -equatorial_j2000[1] * se + equatorial_j2000[2] * ce,
    )


def test_mean_lunar_speeds_are_central_derivatives() -> None:
    """Every lunar node/apsis speed is the derivative of its position slot."""
    flags = eph.FLG_SWIEPH | eph.FLG_SPEED | eph.FLG_TRUEPOS | eph.FLG_NONUT
    current = eph.nod_aps(J2000, eph.MOON, eph.NODBIT_MEAN, flags)
    step = 0.001
    previous = eph.nod_aps(
        J2000 - step,
        eph.MOON,
        eph.NODBIT_MEAN,
        flags & ~eph.FLG_SPEED,
    )
    following = eph.nod_aps(
        J2000 + step,
        eph.MOON,
        eph.NODBIT_MEAN,
        flags & ~eph.FLG_SPEED,
    )
    for point, before, after in zip(current, previous, following):
        expected_longitude_rate = ((after[0] - before[0] + 180.0) % 360.0 - 180.0) / (
            2.0 * step
        )
        assert point[3] == pytest.approx(expected_longitude_rate, abs=2e-10)
        assert point[4] == pytest.approx(
            (after[1] - before[1]) / (2.0 * step), abs=2e-10
        )
        assert point[5] == pytest.approx(
            (after[2] - before[2]) / (2.0 * step), abs=2e-12
        )


@pytest.mark.parametrize("body", [eph.MOON, eph.MERCURY, eph.NEPTUNE])
def test_barycentric_center_uses_compatibility_heliocentric_state(body: int) -> None:
    """The compatibility BARYCTR and HELCTR center selectors are aliases."""
    base = eph.FLG_SWIEPH | eph.FLG_SPEED | eph.FLG_TRUEPOS | eph.FLG_NONUT
    heliocentric = eph.nod_aps(J2000, body, eph.NODBIT_MEAN, base | eph.FLG_HELCTR)
    barycentric = eph.nod_aps(J2000, body, eph.NODBIT_MEAN, base | eph.FLG_BARYCTR)
    for helio_point, bary_point in zip(heliocentric, barycentric):
        assert bary_point == pytest.approx(helio_point, abs=2e-13)


@pytest.mark.parametrize("jd_tt", [2415020.0, J2000, 2460310.5])
def test_output_representation_flags_are_coordinate_identities(jd_tt: float) -> None:
    base = eph.FLG_SWIEPH | eph.FLG_SPEED | eph.FLG_TRUEPOS | eph.FLG_NONUT
    spherical = eph.nod_aps(jd_tt, eph.MARS, eph.NODBIT_MEAN, base)
    cartesian = eph.nod_aps(jd_tt, eph.MARS, eph.NODBIT_MEAN, base | eph.FLG_XYZ)
    radians = eph.nod_aps(jd_tt, eph.MARS, eph.NODBIT_MEAN, base | eph.FLG_RADIANS)
    for sphere, xyz, radian in zip(spherical, cartesian, radians):
        assert xyz == pytest.approx(_spherical_state_to_cartesian(sphere), abs=2e-13)
        assert radian == pytest.approx(
            (
                math.radians(sphere[0]),
                math.radians(sphere[1]),
                sphere[2],
                math.radians(sphere[3]),
                math.radians(sphere[4]),
                sphere[5],
            ),
            abs=2e-13,
        )
        assert all(type(value) is float for value in xyz)
        assert all(type(value) is float for value in radian)


@pytest.mark.parametrize("jd_tt", [2415020.0, J2000, 2460310.5])
def test_equatorial_flag_is_obliquity_rotation(jd_tt: float) -> None:
    base = (
        eph.FLG_SWIEPH
        | eph.FLG_SPEED
        | eph.FLG_HELCTR
        | eph.FLG_TRUEPOS
        | eph.FLG_NONUT
    )
    ecliptic = eph.nod_aps(jd_tt, eph.MARS, eph.NODBIT_MEAN, base | eph.FLG_XYZ)
    equatorial = eph.nod_aps(
        jd_tt,
        eph.MARS,
        eph.NODBIT_MEAN,
        base | eph.FLG_XYZ | eph.FLG_EQUATORIAL,
    )
    eps = vondrak_mean_obliquity_rad(jd_tt)
    ce, se = math.cos(eps), math.sin(eps)
    for ecliptic_state, equatorial_state in zip(ecliptic, equatorial):
        expected = (
            ecliptic_state[0],
            ecliptic_state[1] * ce - ecliptic_state[2] * se,
            ecliptic_state[1] * se + ecliptic_state[2] * ce,
            ecliptic_state[3],
            ecliptic_state[4] * ce - ecliptic_state[5] * se,
            ecliptic_state[4] * se + ecliptic_state[5] * ce,
        )
        assert equatorial_state == pytest.approx(expected, abs=2e-12)


@pytest.mark.parametrize(
    ("jd_tt", "body", "center_flag"),
    [
        (2415020.0, eph.MARS, eph.FLG_HELCTR),
        (J2000, eph.MOON, 0),
        (2460310.5, eph.NEPTUNE, eph.FLG_HELCTR),
    ],
)
def test_j2000_speed_is_derivative_of_time_dependent_frame_rotation(
    jd_tt: float, body: int, center_flag: int
) -> None:
    base = (
        eph.FLG_SWIEPH | eph.FLG_SPEED | eph.FLG_TRUEPOS | eph.FLG_NONUT | center_flag
    )
    of_date = eph.nod_aps(jd_tt, body, eph.NODBIT_MEAN, base)
    j2000 = eph.nod_aps(jd_tt, body, eph.NODBIT_MEAN, base | eph.FLG_J2000)
    frame_step = 0.001
    for date_point, j2000_point in zip(of_date, j2000):
        date_state = _spherical_state_to_cartesian(date_point)
        expected_position = _ecliptic_date_to_j2000(date_state[:3], jd_tt)
        rotated_velocity = _ecliptic_date_to_j2000(date_state[3:], jd_tt)
        before = _ecliptic_date_to_j2000(date_state[:3], jd_tt - frame_step)
        after = _ecliptic_date_to_j2000(date_state[:3], jd_tt + frame_step)
        expected_velocity = tuple(
            rotated_velocity[index]
            + (after[index] - before[index]) / (2.0 * frame_step)
            for index in range(3)
        )
        actual = _spherical_state_to_cartesian(j2000_point)
        assert actual[:3] == pytest.approx(expected_position, abs=2e-12)
        assert actual[3:] == pytest.approx(expected_velocity, abs=2e-10)


@pytest.mark.parametrize("jd_tt", [2415020.0, J2000, 2460310.5])
def test_topocentric_true_state_subtracts_iau_observer_state(jd_tt: float) -> None:
    eph.set_topo(12.5, 41.9, 100.0)
    flags = eph.FLG_SWIEPH | eph.FLG_SPEED | eph.FLG_TRUEPOS | eph.FLG_NONUT
    geocentric = eph.nod_aps(jd_tt, eph.MARS, eph.NODBIT_MEAN, flags)
    topocentric = eph.nod_aps(jd_tt, eph.MARS, eph.NODBIT_MEAN, flags | eph.FLG_TOPOCTR)

    ts = get_timescale()
    epoch = ts.tt_jd(jd_tt)
    earth = get_planets()["earth"]
    observer = get_topo()
    assert observer is not None
    topocenter = (earth + observer).at(epoch)
    geocenter = earth.at(epoch)
    position_icrs = tuple(
        float(topocenter.position.au[index] - geocenter.position.au[index])
        for index in range(3)
    )
    velocity_icrs = tuple(
        float(topocenter.velocity.au_per_d[index] - geocenter.velocity.au_per_d[index])
        for index in range(3)
    )
    expected_position = _icrs_to_ecliptic(position_icrs, jd_tt, nutation=False)
    expected_velocity = _icrs_to_ecliptic(velocity_icrs, jd_tt, nutation=False)
    frame_step = 0.001
    before = _icrs_to_ecliptic(position_icrs, jd_tt - frame_step, nutation=False)
    after = _icrs_to_ecliptic(position_icrs, jd_tt + frame_step, nutation=False)
    expected_velocity = tuple(
        expected_velocity[index] + (after[index] - before[index]) / (2.0 * frame_step)
        for index in range(3)
    )

    for geo_point, topo_point in zip(geocentric, topocentric):
        geo_state = _spherical_state_to_cartesian(geo_point)
        topo_state = _spherical_state_to_cartesian(topo_point)
        offset = tuple(geo_state[index] - topo_state[index] for index in range(6))
        assert offset[:3] == pytest.approx(expected_position, abs=2e-13)
        assert offset[3:] == pytest.approx(expected_velocity, abs=2e-10)


def test_topocentric_diurnal_aberration_uses_rotation_velocity() -> None:
    eph.set_topo(12.5, 41.9, 100.0)
    flags = eph.FLG_SWIEPH | eph.FLG_NOGDEFL
    geocentric = eph.nod_aps(J2000, eph.MARS, eph.NODBIT_MEAN, flags)
    topocentric = eph.nod_aps(J2000, eph.MARS, eph.NODBIT_MEAN, flags | eph.FLG_TOPOCTR)

    ts = get_timescale()
    epoch = ts.tt_jd(J2000)
    earth = get_planets()["earth"]
    observer = get_topo()
    assert observer is not None
    topocenter = (earth + observer).at(epoch)
    geocenter = earth.at(epoch)
    offset = _icrs_to_ecliptic(
        tuple(
            float(topocenter.position.au[index] - geocenter.position.au[index])
            for index in range(3)
        ),
        J2000,
        nutation=True,
    )
    relative_velocity = _icrs_to_ecliptic(
        tuple(
            float(
                topocenter.velocity.au_per_d[index] - geocenter.velocity.au_per_d[index]
            )
            for index in range(3)
        ),
        J2000,
        nutation=True,
    )

    for geo_point, topo_point in zip(geocentric, topocentric):
        geo_position = _spherical_state_to_cartesian(geo_point)[:3]
        expected = apply_aberration_to_position(
            tuple(geo_position[index] - offset[index] for index in range(3)),
            relative_velocity,
        )
        actual = _spherical_state_to_cartesian(topo_point)[:3]
        assert actual == pytest.approx(expected, abs=2e-13)


def test_topocentric_without_speed_returns_zero_speed_slots() -> None:
    eph.set_topo(12.5, 41.9, 100.0)
    points = eph.nod_aps(
        J2000,
        eph.MARS,
        eph.NODBIT_MEAN,
        eph.FLG_SWIEPH | eph.FLG_TOPOCTR,
    )
    assert all(point[3:] == (0.0, 0.0, 0.0) for point in points)


def test_topocentric_preserves_undefined_zero_node_sentinels() -> None:
    eph.set_topo(12.5, 41.9, 100.0)
    ascending, descending, _, _ = eph.nod_aps(
        J2000,
        eph.SUN,
        eph.NODBIT_MEAN,
        eph.FLG_SWIEPH | eph.FLG_SPEED | eph.FLG_TOPOCTR,
    )
    assert ascending == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert descending == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def test_apparent_place_suppression_flags_compose_without_recorded_fixtures() -> None:
    base = eph.FLG_SWIEPH | eph.FLG_SPEED
    geometric = eph.nod_aps(J2000, eph.MARS, eph.NODBIT_MEAN, base | eph.FLG_TRUEPOS)
    suppressed = eph.nod_aps(
        J2000,
        eph.MARS,
        eph.NODBIT_MEAN,
        base | eph.FLG_NOABERR | eph.FLG_NOGDEFL,
    )
    for geometric_point, suppressed_point in zip(geometric, suppressed):
        assert suppressed_point == pytest.approx(geometric_point, abs=2e-12)
