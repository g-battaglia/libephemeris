"""Public values share the native engine's scalar and finite table definitions."""

import inspect
import json
from pathlib import Path

import libephemeris as ephem
from libephemeris import _core
from libephemeris import constants
import pytest


def test_numeric_constant_types_and_canonical_aliases():
    assert type(constants.SUN) is int
    assert type(constants.J2000) is float
    assert type(constants.SE_FNAME_DE431) is str
    assert constants.SE_FNAME_DE431 == constants.FNAME_DE431
    assert constants.TRANSPLUTO == constants.ISIS
    assert constants.WHITE_MOON == constants.SELENA
    assert constants.AUNIT_TO_KM == 149597870.7


def test_finite_map_layouts_and_native_lookup_consistency():
    assert type(_core.INNER_PLANETS) is set
    assert _core.INNER_PLANETS == {constants.MERCURY, constants.VENUS}
    assert _core.EXPERIENCE_FACTORS == {1: 1.3, 2: 1.0, 3: 0.85, 4: 0.7}
    assert type(constants.REQUIRED_SPK_BODIES) is set
    for body, pair in constants.SPK_BODY_NAME_MAP.items():
        horizons, naif = pair
        assert type(body) is int and type(horizons) is str and type(naif) is int
        assert constants.get_horizons_id(body) == horizons
        assert constants.get_naif_id_from_ipl(body) == naif
        assert constants.get_spk_body_info_from_map(body) == pair
        assert constants.is_spk_auto_download_blocked(body) == (
            body in constants.SPK_AUTO_DOWNLOAD_BLOCKED
        )
    assert all(
        body in constants.SPK_BODY_NAME_MAP for body in constants.REQUIRED_SPK_BODIES
    )
    assert _core.MOON_NAIF_MAP[constants.MOON_IO] == constants.NAIF_IO
    assert _core.MOON_NAMES[constants.MOON_CHARON] == "Charon"
    assert _core.MOON_PARENT_MAP[constants.MOON_CHARON] == constants.PLUTO
    assert _core.MAJOR_ASTEROID_SPK_INFO[constants.CERES] == (
        1,
        "Ceres;",
        constants.NAIF_CERES,
        "Ceres",
    )


@pytest.mark.parametrize("body", [-2147483648, -1, 0, 2147483647])
def test_unknown_lookup_keys(body):
    assert constants.get_horizons_id(body) is None
    assert constants.get_naif_id_from_ipl(body) is None
    assert constants.get_spk_body_info_from_map(body) is None
    assert constants.is_spk_auto_download_blocked(body) is False


@pytest.mark.parametrize(
    "name",
    [
        "get_horizons_id",
        "get_naif_id_from_ipl",
        "get_spk_body_info_from_map",
        "is_spk_auto_download_blocked",
    ],
)
def test_lookup_signatures_and_standard_integer_conversion(name):
    function = getattr(constants, name)
    signature = inspect.signature(function)
    assert list(signature.parameters) == ["ipl"]
    assert signature.parameters["ipl"].default is inspect.Parameter.empty
    with pytest.raises(TypeError):
        function("15")
