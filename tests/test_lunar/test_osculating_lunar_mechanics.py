# SPDX-License-Identifier: Apache-2.0
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
