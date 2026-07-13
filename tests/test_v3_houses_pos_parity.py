"""Reference-free invariants for high-latitude ``house_pos`` systems."""

from __future__ import annotations

import importlib

import pytest


H = importlib.import_module("libephemeris.houses")

CASES = [
    (123.0, 66.0, 23.4393, 47.3, 2.5),
    (271.0, 71.0, 23.4393, 251.9, -4.5),
    (34.0, -63.0, 23.4393, 133.6, 4.5),
    (317.0, 68.5, 23.4368, 299.3, -3.0),
    (200.0, -72.0, 23.4393, 88.1, 0.0),
    (95.0, 74.0, 23.437, 331.7, 5.0),
]
SYSTEMS = "EAWVTHFYIiJSU"


def _house_delta(a: float, b: float) -> float:
    delta = abs(a - b)
    return min(delta, 12.0 - delta)


@pytest.mark.parametrize("hsys", SYSTEMS)
@pytest.mark.parametrize("armc,geolat,eps,lon,lat_body", CASES)
def test_high_latitude_results_are_finite_bounded_and_periodic(
    hsys: str,
    armc: float,
    geolat: float,
    eps: float,
    lon: float,
    lat_body: float,
) -> None:
    got = H.house_pos(armc, geolat, eps, [lon, lat_body], hsys.encode())
    wrapped = H.house_pos(
        armc + 360.0, geolat, eps, [lon + 360.0, lat_body], hsys.encode()
    )
    assert 1.0 <= got < 13.0
    assert _house_delta(got, wrapped) < 1e-10


@pytest.mark.parametrize("hsys", ["E", "A", "W", "V"])
def test_equal_systems_advance_one_house_per_thirty_degrees(hsys: str) -> None:
    first = H.house_pos(330.0, 75.0, 23.4393, [300.0, 0.0], hsys.encode())
    second = H.house_pos(330.0, 75.0, 23.4393, [330.0, 0.0], hsys.encode())
    assert ((second - first) % 12.0) == pytest.approx(1.0, abs=1e-9)


def test_topocentric_uses_a_distinct_position_circle() -> None:
    armc, geolat, eps, lon, lat_body = 271.0, 71.0, 23.4393, 251.9, -4.5
    topocentric = H.house_pos(armc, geolat, eps, [lon, lat_body], b"T")
    placidus = H.house_pos(armc, geolat, eps, [lon, lat_body], b"P")
    assert _house_delta(topocentric, placidus) > 0.01


def test_horizon_system_depends_on_ecliptic_latitude() -> None:
    armc, geolat, eps, lon = 123.0, 66.0, 23.4393, 47.3
    north = H.house_pos(armc, geolat, eps, [lon, 5.0], b"H")
    equator = H.house_pos(armc, geolat, eps, [lon, 0.0], b"H")
    south = H.house_pos(armc, geolat, eps, [lon, -5.0], b"H")
    assert north != equator
    assert south != equator


@pytest.mark.parametrize("hsys", ["I", "i", "Y"])
def test_sunshine_and_apc_depend_on_latitude(hsys: str) -> None:
    top = H.house_pos(271.0, 71.0, 23.4393, [251.9, 6.0], hsys.encode())
    bottom = H.house_pos(271.0, 71.0, 23.4393, [251.9, -6.0], hsys.encode())
    assert top != bottom
