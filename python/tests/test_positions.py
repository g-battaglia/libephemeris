"""Typed position forwarding and concurrent immutable snapshots."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
import inspect
from pathlib import Path

import libephemeris as ephem
import pytest


@pytest.fixture(autouse=True)
def explicit_lephemeris():
    ephem.set_calc_mode("leb")
    ephem.set_precision_tier("base")
    ephem.set_leb_file(str(CORE_LEB))
    ephem.set_iers_delta_t_enabled(False)
    ephem.set_delta_t_userdef(64.0 / 86400.0)


def test_typed_layout_and_signatures():
    for name, epoch_parameter in (("calc_ut", "tjdut"), ("calc", "tjdet")):
        function = getattr(ephem, name)
        signature = inspect.signature(function)
        assert list(signature.parameters) == [epoch_parameter, "planet", "flags"]
        assert signature.parameters["flags"].default == 258
        values, flags = function(**{epoch_parameter: 2451545.0, "planet": 0})
        assert len(values) == 6 and all(type(value) is float for value in values)
        assert type(flags) is int
        with pytest.raises(TypeError):
            function(object(), 0)


def test_eight_threads_equal_sequential_calls():
    requests = [(2451545.0 + day, body, 258) for day in range(4) for body in range(10)]
    expected = [ephem.calc_ut(*request) for request in requests]
    with ThreadPoolExecutor(max_workers=8) as executor:
        actual = list(executor.map(lambda request: ephem.calc_ut(*request), requests))
    assert actual == expected


def test_invalid_configuration_is_value_error():
    with pytest.raises(ValueError):
        ephem.set_calc_mode("invalid")
    with pytest.raises(ValueError):
        ephem.set_precision_tier("invalid")


@pytest.mark.parametrize("api", ["calc", "calc_ut"])
@pytest.mark.parametrize("epoch", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_epochs_are_refused_at_the_public_boundary(api, epoch):
    with pytest.raises(ValueError):
        getattr(ephem, api)(epoch, ephem.MARS)
