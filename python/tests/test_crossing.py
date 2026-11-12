"""Physical roots, direction conversion and native crossing signatures."""
from data_paths import CORE_LEB

from pathlib import Path
import inspect

import libephemeris as eph
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_iers_delta_t_enabled(False)


@pytest.mark.parametrize("function,position,body", [
    ("solcross_ut", "calc_ut", 0), ("solcross", "calc", 0),
    ("mooncross_ut", "calc_ut", 1), ("mooncross", "calc", 1),
])
def test_longitude_roots_and_directions(function, position, body):
    search, calculate = getattr(eph, function), getattr(eph, position)
    for backwards in (False, True):
        root = search(90.0, 2451545.0, 2, backwards)
        assert root < 2451545.0 if backwards else root >= 2451545.0
        values, _ = calculate(root, body, 258)
        residual = ((values[0] - 90.0 + 180.0) % 360.0) - 180.0
        assert abs(residual) * 3600 <= 0.001
    signature = inspect.signature(search)
    assert signature.parameters["flags"].default == 2
    assert signature.parameters["backwards"].default is False


def test_node_and_station_layouts():
    time, longitude, latitude = eph.mooncross_node_ut(2451545.0)
    assert time >= 2451545.0
    assert 0 <= longitude < 360
    assert abs(latitude) * 3600 <= 0.001
    time, kind = eph.find_station_ut(2, 2451545.0)
    assert time > 2451545.0 and kind in ("SR", "SD")
    start, end = eph.next_retrograde_ut(2, 2451545.0)
    assert 2451545.0 < start < end


def test_heliocentric_root_uses_native_position_pipeline():
    time = eph.helio_cross_ut(2, 90.0, 2451545.0)
    values, _ = eph.calc_ut(time, 2, 2 | 8 | 256)
    residual = ((values[0] - 90.0 + 180.0) % 360.0) - 180.0
    assert time >= 2451545.0 and abs(residual) * 3600 <= 0.001


@pytest.mark.parametrize("direction,backwards", [("backward", True), ("forward", False), (1, True), (0, False)])
def test_documented_direction_types(direction, backwards):
    assert eph.solcross_ut(0.0, 2451545.0, 2, direction) == eph.solcross_ut(0.0, 2451545.0, 2, backwards)


def test_invalid_direction_and_station_body():
    with pytest.raises(TypeError):
        eph.solcross_ut(0.0, 2451545.0, backwards="invalid")
    with pytest.raises(TypeError):
        eph.solcross_ut(0.0, 2451545.0, backwards=None)
    with pytest.raises(ValueError):
        eph.find_station_ut(0, 2451545.0)
