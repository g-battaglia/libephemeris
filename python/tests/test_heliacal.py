"""Native visibility layouts, physical predicates and event input domains."""
from data_paths import CORE_LEB, LEB_HASHES

import hashlib
import inspect
import math
from pathlib import Path

import libephemeris as eph
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_iers_delta_t_enabled(False)


def test_body_classification_uses_native_integer_domain():
    assert eph.is_inner_planet(2) and eph.is_inner_planet(3)
    assert not any(eph.is_inner_planet(body) for body in (-1, 0, 1, 4, 14, 1000000))
    assert not eph.is_fixed_star(999999)
    assert eph.is_fixed_star(1000000) and eph.is_fixed_star(1000001)


def test_visibility_result_layout_and_reserved_slots():
    flag, values = eph.vis_limit_mag(2451545.0, (12.5, 41.9, 30.0), (), (), "Venus")
    assert isinstance(flag, float)
    assert isinstance(values, tuple) and len(values) == 10
    assert all(math.isfinite(value) for value in values)
    assert values[8:] == (0.0, 0.0)
    if flag == -2.0:
        assert values == (-100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    else:
        assert -90 <= values[1] <= 90
        assert 0 <= values[2] < 360


def test_signatures_and_typed_inputs():
    for name in ("heliacal_ut", "heliacal_pheno_ut", "vis_limit_mag"):
        signature = inspect.signature(getattr(eph, name))
        assert signature.parameters["flags"].default == 2
        assert "objname" in signature.parameters
    with pytest.raises(TypeError):
        eph.vis_limit_mag(2451545.0, (), (), (), object())
    with pytest.raises(eph.Error):
        eph.heliacal_ut(2451545.0, (), (), (), "Venus", 99)


def test_native_pheno_keeps_fifty_physical_slots():
    from data_paths import DATA_ROOT
    kernel = DATA_ROOT / "de440.bsp"
    if not kernel.is_file():
        pytest.skip("optional DE440 kernel unavailable")
    eph.set_calc_mode("skyfield")
    eph.set_leb_file(None)
    eph.set_ephe_path(str(kernel.parent))
    eph.set_ephemeris_file(str(kernel))
    values = eph.heliacal_pheno_ut(2451545.0, (12.5, 41.9, 30.0), (), (), "Venus", 1)
    assert isinstance(values, tuple) and len(values) == 50
    assert all(math.isfinite(value) for value in values)
    assert values[28:] == (0.0,) * 22
    assert 0 <= values[27] <= 100


@pytest.mark.parametrize("mode", ["leb", "skyfield"])
def test_mercury_setting_returns_a_positive_public_visibility_interval(mode):
    from data_paths import DATA_ROOT, LEB_HASHES

    core = DATA_ROOT / "leb/medium_core.leb2"
    kernel = DATA_ROOT / "de440.bsp"
    if not core.is_file() or not kernel.is_file():
        pytest.skip("optional medium source inputs unavailable")
    with core.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == LEB_HASHES["medium"]["core"]
    with kernel.open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"
    eph.set_network_policy("sealed")
    eph.set_auto_spk_download(False)
    eph.set_precision_tier("medium")
    eph.set_calc_mode(mode)
    eph.set_leb_file(str(core) if mode == "leb" else None)
    eph.set_ephe_path(str(DATA_ROOT))
    eph.set_ephemeris_file(str(kernel))
    eph.set_delta_t_userdef(None)
    arguments = ((0.0, 30.0, 0.0), (1013.25, 15.0, 50.0, 0.25),
                 (25.0, 1.0, 1.0, 1.0, 0.0, 0.0), "Mercury")
    first, optimum, last = eph.heliacal_ut(2460310.5, *arguments, 2)
    assert 2460310.5 <= first < optimum < last

    def margin(epoch):
        status, values = eph.vis_limit_mag(epoch, *arguments)
        assert all(math.isfinite(value) for value in (status, *values))
        return -99.0 if status < 0 or values[3] >= 0 else values[0] - values[7]

    assert margin(optimum) > 0.0
    second = 1.0 / 86400.0
    assert margin(first - second) <= 0.0 < margin(first + second)
    assert margin(last + second) <= 0.0 < margin(last - second)
