# SPDX-License-Identifier: AGPL-3.0-only
"""Tests for the JPL-state-vector lunar node and osculating apsides."""

from __future__ import annotations

import math

import pytest
from skyfield.framelib import ecliptic_frame

from libephemeris.cache import get_cached_time_tt
from libephemeris.lunar import (
    GM_EARTH_MOON_AU3_DAY2,
    calc_osculating_perigee,
    calc_true_lilith,
    calc_true_lilith_orbital_elements,
    calc_true_lunar_node,
)
from libephemeris.state import get_planets

_GM_EARTH_KM3_S2 = 398600.435436
_EARTH_MOON_MASS_RATIO = 81.3005691
_AU_KM = 149597870.7
_SECONDS_PER_DAY = 86400.0


def _cross(
    left: tuple[float, float, float], right: tuple[float, float, float]
) -> tuple[float, float, float]:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(sum(component * component for component in vector))


def _moon_state(
    jd_tt: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    planets = get_planets()
    time = get_cached_time_tt(jd_tt)
    position = (planets["moon"] - planets["earth"]).at(time)
    xyz, velocity = position.frame_xyz_and_velocity(ecliptic_frame)
    return (
        tuple(float(value) for value in xyz.au),
        tuple(float(value) for value in velocity.au_per_d),
    )


def _eccentricity_vector(
    position: tuple[float, float, float],
    velocity: tuple[float, float, float],
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    angular_momentum = _cross(position, velocity)
    velocity_cross_momentum = _cross(velocity, angular_momentum)
    radius = _norm(position)
    eccentricity = tuple(
        velocity_cross_momentum[index] / GM_EARTH_MOON_AU3_DAY2
        - position[index] / radius
        for index in range(3)
    )
    return eccentricity, angular_momentum  # type: ignore[return-value]


def _spherical(vector: tuple[float, float, float]) -> tuple[float, float]:
    magnitude = _norm(vector)
    return (
        float(math.degrees(math.atan2(vector[1], vector[0])) % 360.0),
        float(math.degrees(math.asin(vector[2] / magnitude))),
    )


def test_combined_earth_moon_gravitational_parameter() -> None:
    moon_gm = _GM_EARTH_KM3_S2 / _EARTH_MOON_MASS_RATIO
    expected = (_GM_EARTH_KM3_S2 + moon_gm) / _AU_KM**3 * _SECONDS_PER_DAY**2
    assert GM_EARTH_MOON_AU3_DAY2 == pytest.approx(expected, rel=1e-15)


@pytest.mark.parametrize("jd_tt", [2451545.0, 2458849.5, 2460676.5])
def test_osculating_apsides_follow_the_eccentricity_vector(jd_tt: float) -> None:
    position, velocity = _moon_state(jd_tt)
    eccentricity, angular_momentum = _eccentricity_vector(position, velocity)
    eccentricity_magnitude = _norm(eccentricity)
    semilatus_rectum = _norm(angular_momentum) ** 2 / GM_EARTH_MOON_AU3_DAY2

    expected_perigee = (
        *_spherical(eccentricity),
        semilatus_rectum / (1.0 + eccentricity_magnitude),
    )
    expected_apogee = (
        *_spherical(tuple(-component for component in eccentricity)),
        semilatus_rectum / (1.0 - eccentricity_magnitude),
    )

    perigee = calc_osculating_perigee(jd_tt)
    apogee = calc_true_lilith(jd_tt)
    assert perigee == pytest.approx(expected_perigee, abs=2e-14)
    assert apogee == pytest.approx(expected_apogee, abs=2e-14)
    assert all(type(value) is float for value in perigee + apogee)


@pytest.mark.parametrize("jd_tt", [2451545.0, 2458849.5, 2460676.5])
def test_osculating_apogee_and_perigee_are_opposite(jd_tt: float) -> None:
    apogee = calc_true_lilith(jd_tt)
    perigee = calc_osculating_perigee(jd_tt)
    separation = (perigee[0] - apogee[0]) % 360.0

    assert separation == pytest.approx(180.0, abs=2e-12)
    assert perigee[1] == pytest.approx(-apogee[1], abs=2e-12)
    assert 0.002 < perigee[2] < apogee[2] < 0.004
    assert calc_true_lilith_orbital_elements(jd_tt) == apogee


@pytest.mark.parametrize(
    "jd_tt",
    [
        2413687.5,  # 1896
        2426405.5,  # 1931
        2435130.5,  # 1955
        2447515.5,  # 1988
        2459034.5,  # 2020
        2478361.5,  # 2073
    ],
)
def test_osculating_apogee_speed_is_self_consistent(jd_tt: float) -> None:
    """The FLG_SPEED longitude speed reported for the osculating (true) apogee
    (OscuApog, id 13) must be the genuine time-derivative of the reported
    longitude curve. A coarse central-difference half-step chords across the
    apogee's fast short-period structure and biased this speed by ~0.2-1.2"/day;
    the finer half-step keeps it self-consistent to well under 0.03"/day.
    """
    import libephemeris as ephem
    from libephemeris.constants import OSCU_APOG, FLG_SPEED

    def _longitude(jd: float) -> float:
        return calc_true_lilith(jd)[0]

    def _unwrap(delta: float) -> float:
        if delta > 180.0:
            return delta - 360.0
        if delta < -180.0:
            return delta + 360.0
        return delta

    # Independent 4th-order central derivative of the reported longitude curve.
    h = 0.002
    d1 = _unwrap(_longitude(jd_tt + h) - _longitude(jd_tt - h))
    d2 = _unwrap(_longitude(jd_tt + 2.0 * h) - _longitude(jd_tt - 2.0 * h))
    true_speed = (8.0 * d1 - d2) / (12.0 * h)

    position, _ = ephem.calc(jd_tt, OSCU_APOG, FLG_SPEED)
    reported_speed = position[3]

    auto_inconsistency_arcsec = abs(reported_speed - true_speed) * 3600.0
    assert auto_inconsistency_arcsec < 0.03, (
        f"OscuApog longitude speed {reported_speed:.6f} deg/day is "
        f"self-inconsistent with its own longitude derivative "
        f'{true_speed:.6f} deg/day by {auto_inconsistency_arcsec:.4f} "/day'
    )


@pytest.mark.parametrize("jd_tt", [2451545.0, 2458849.5, 2460676.5])
def test_true_node_follows_the_orbital_plane_intersection(jd_tt: float) -> None:
    position, velocity = _moon_state(jd_tt)
    eccentricity, angular_momentum = _eccentricity_vector(position, velocity)
    h_x, h_y, h_z = angular_momentum
    expected_longitude = math.degrees(math.atan2(h_x, -h_y)) % 360.0

    node_magnitude = math.hypot(h_x, h_y)
    eccentricity_magnitude = _norm(eccentricity)
    cos_argument_of_perigee = (-h_y * eccentricity[0] + h_x * eccentricity[1]) / (
        node_magnitude * eccentricity_magnitude
    )
    argument_of_perigee = math.acos(max(-1.0, min(1.0, cos_argument_of_perigee)))
    if eccentricity[2] < 0.0:
        argument_of_perigee = 2.0 * math.pi - argument_of_perigee
    semilatus_rectum = (h_x * h_x + h_y * h_y + h_z * h_z) / GM_EARTH_MOON_AU3_DAY2
    expected_distance = semilatus_rectum / (
        1.0 + eccentricity_magnitude * math.cos(2.0 * math.pi - argument_of_perigee)
    )

    longitude, latitude, distance = calc_true_lunar_node(jd_tt)
    assert longitude == pytest.approx(expected_longitude, abs=2e-12)
    assert latitude == 0.0
    assert distance == pytest.approx(expected_distance, abs=2e-14)
    assert all(type(value) is float for value in (longitude, latitude, distance))


def test_interpolated_apsides_speed_is_self_consistent():
    """INTP_APOG/INTP_PERG FLG_SPEED equals the derivative of the reported
    positions: the half-day stencil misrepresented the short-period
    Delaunay/residual structure by up to ~10"/day at fast-swing phases."""
    import libephemeris as le

    cases = [(2079, 9, 17, 2.4, 22), (1940, 6, 15, 12.0, 22), (2001, 1, 21, 2.4, 21)]
    f = le.FLG_SWIEPH | le.FLG_SPEED
    for y, m, d, h, body in cases:
        jd = le.julday(y, m, d, h)
        rep = le.calc_ut(jd, body, f)[0]
        hh = 0.02
        p_m = le.calc_ut(jd - hh, body, f)[0]
        p_p = le.calc_ut(jd + hh, body, f)[0]
        dl = (p_p[0] - p_m[0] + 180.0) % 360.0 - 180.0
        assert abs(rep[3] - dl / (2 * hh)) * 3600.0 < 0.3
