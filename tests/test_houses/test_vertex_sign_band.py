"""Baseline-free projective tests for the auxiliary ASCMC geometry."""

from __future__ import annotations

import math

import pytest

import libephemeris as le
from libephemeris.houses import _auxiliary_ascmc, _calc_vertex

EPS = 23.4393


def _vertex(lat: float, hsys: str) -> float:
    return le.houses_armc(0.0, lat, EPS, ord(hsys))[1][3]


def _angular_difference(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def test_placidus_one_sided_limits_are_antipodal():
    assert _angular_difference(_vertex(-1e-11, "P"), _vertex(1e-11, "P")) > 179.0


def test_exact_zero_selects_placidus_positive_limit():
    assert _angular_difference(_vertex(0.0, "P"), _vertex(1e-11, "P")) < 1e-6


def test_exact_zero_selects_horizontal_positive_limit():
    assert _angular_difference(_vertex(0.0, "H"), _vertex(1e-11, "H")) < 1e-6


@pytest.mark.parametrize(
    ("armc", "latitude", "expected"),
    [(90.0, EPS, 90.0), (270.0, -EPS, 270.0)],
)
def test_ordinary_plane_coincidences_use_right_hand_limit(
    armc: float, latitude: float, expected: float
):
    assert _calc_vertex(armc, EPS, latitude) == expected


def test_zero_west_coordinate_uses_analytic_armc_derivative():
    armc = 30.0
    eps = 45.0
    sin_theta = math.sin(math.radians(armc))
    cos_theta = math.cos(math.radians(armc))
    sin_eps = math.sin(math.radians(eps))
    cos_eps = math.cos(math.radians(eps))
    latitude = math.degrees(math.atan2(sin_eps * sin_theta, cos_eps))
    sin_latitude = math.sin(math.radians(latitude))
    cos_latitude = math.cos(math.radians(latitude))

    branch = sin_eps * cos_latitude * sin_theta - cos_eps * sin_latitude
    derivative = sin_eps * cos_latitude * cos_theta
    x = cos_eps * sin_latitude * sin_theta - sin_eps * cos_latitude
    y = -sin_latitude * cos_theta
    assert branch == 0.0
    assert derivative > 0.0

    expected = math.degrees(math.atan2(-y, -x)) % 360.0
    assert _calc_vertex(armc, eps, latitude) == expected


@pytest.mark.parametrize("armc", [0.0, 47.0, 90.0, 271.0, 360.0])
def test_zero_obliquity_vertex_is_west_point(armc: float):
    assert _calc_vertex(armc, 0.0, 38.0) == (armc - 90.0) % 360.0


@pytest.mark.parametrize(
    ("latitude", "armc", "expected"),
    [
        (0.0, 0.0, 180.0),
        (0.0, 90.0, 0.0),
        (0.0, 180.0, 180.0),
        (0.0, 270.0, 180.0),
        (90.0, 0.0, 270.0),
        (90.0, 90.0, 0.0),
        (90.0, 180.0, 90.0),
        (90.0, 270.0, 180.0),
        (-90.0, 0.0, 270.0),
        (-90.0, 90.0, 0.0),
        (-90.0, 180.0, 90.0),
        (-90.0, 270.0, 180.0),
    ],
)
def test_right_angle_obliquity_exact_cells(
    latitude: float, armc: float, expected: float
):
    assert _calc_vertex(armc, 90.0, latitude) == expected


@pytest.mark.parametrize("latitude", [-90.0, 90.0])
def test_plane_defined_slots_at_geographic_poles(latitude: float):
    _, east, koch, munkasey, polar = _auxiliary_ascmc(37.0, latitude, EPS, "P")
    assert polar == (180.0 if latitude > 0.0 else 0.0)
    assert koch == (0.0 if latitude > 0.0 else 180.0)
    assert munkasey == east


def test_plane_defined_slot_relationships():
    armc = 37.0
    latitude = 28.0
    _, east, koch, munkasey, polar = _auxiliary_ascmc(armc, latitude, EPS, "P")

    east_ra = (
        math.degrees(
            math.atan2(
                math.cos(math.radians(EPS)) * math.sin(math.radians(east)),
                math.cos(math.radians(east)),
            )
        )
        % 360.0
    )
    assert _angular_difference(east_ra, armc + 90.0) < 1e-12
    assert (polar + 180.0) % 360.0 == koch

    munkasey_pole = 90.0 - latitude
    alpha = math.radians(armc + 90.0)
    pole = math.radians(munkasey_pole)
    expected_munkasey = (
        math.degrees(
            math.atan2(
                math.sin(alpha) * math.cos(pole),
                math.cos(math.radians(EPS)) * math.cos(alpha) * math.cos(pole)
                - math.sin(math.radians(EPS)) * math.sin(pole),
            )
        )
        % 360.0
    )
    assert munkasey == expected_munkasey


def test_adjacent_binary64_latitudes_are_not_zero_conventions():
    north = math.nextafter(0.0, math.inf)
    south = math.nextafter(0.0, -math.inf)
    placidus_zero = _auxiliary_ascmc(100.0, 0.0, EPS, "P")
    horizontal_zero = _auxiliary_ascmc(100.0, 0.0, EPS, "H")
    north_values = _auxiliary_ascmc(100.0, north, EPS, "H")
    south_values = _auxiliary_ascmc(100.0, south, EPS, "P")

    assert placidus_zero[3] == 180.0
    assert horizontal_zero[3] == 0.0
    assert north_values[3] == 180.0
    assert south_values[3] == 0.0
    assert _angular_difference(north_values[3], south_values[3]) == 180.0


def test_auxiliary_outputs_are_native_normalized_and_periodic():
    values = _auxiliary_ascmc(17.0, -42.0, EPS, "P")
    assert values == _auxiliary_ascmc(17.0 + 5.0 * 360.0, -42.0, EPS, "P")
    assert all(type(value) is float for value in values)
    assert all(0.0 <= value < 360.0 for value in values)


def test_regular_vertex_ray_belongs_to_both_planes_and_is_western():
    armc = math.radians(123.0)
    latitude = math.radians(41.0)
    eps = math.radians(EPS)
    east = (-math.sin(armc), math.cos(armc), 0.0)
    north = (
        -math.sin(latitude) * math.cos(armc),
        -math.sin(latitude) * math.sin(armc),
        math.cos(latitude),
    )
    ecliptic_pole = (0.0, -math.sin(eps), math.cos(eps))
    intersection = (
        ecliptic_pole[1] * north[2] - ecliptic_pole[2] * north[1],
        ecliptic_pole[2] * north[0] - ecliptic_pole[0] * north[2],
        ecliptic_pole[0] * north[1] - ecliptic_pole[1] * north[0],
    )
    west_coordinate = sum(a * b for a, b in zip(intersection, east))
    ray = (
        tuple(-value for value in intersection)
        if west_coordinate > 0.0
        else intersection
    )

    assert sum(a * b for a, b in zip(ray, north)) == pytest.approx(0.0, abs=1e-16)
    assert sum(a * b for a, b in zip(ray, ecliptic_pole)) == pytest.approx(
        0.0, abs=1e-16
    )
    assert sum(a * b for a, b in zip(ray, east)) < 0.0

    longitude = (
        math.degrees(
            math.atan2(
                ray[1] * math.cos(eps) + ray[2] * math.sin(eps),
                ray[0],
            )
        )
        % 360.0
    )
    assert _calc_vertex(math.degrees(armc), EPS, math.degrees(latitude)) == longitude


@pytest.mark.unit
class TestHousePosIntHsysObjcoordForm:
    """The objcoord-first house_pos form honors an int house code (character
    code, same convention as the 6-arg form) instead of silently falling
    back to Placidus."""

    def test_int_matches_bytes_and_sixarg(self):
        from libephemeris.houses import house_pos

        i = house_pos(0, 23.4, 23.4393, (40.0, 0.0), ord("S"))
        b = house_pos(0, 23.4, 23.4393, (40.0, 0.0), b"S")
        six = house_pos(0, 23.4, 23.4393, ord("S"), 40.0, 0.0)
        assert i == b == six
        # And it is NOT the Placidus value.
        p = house_pos(0, 23.4, 23.4393, (40.0, 0.0), b"P")
        assert i != p

    def test_omitted_hsys_defaults_placidus(self):
        from libephemeris.houses import house_pos

        d = house_pos(0, 23.4, 23.4393, (40.0, 0.0))
        p = house_pos(0, 23.4, 23.4393, (40.0, 0.0), b"P")
        assert d == p
