"""Satellite registry ownership, physical outputs and public calc integration."""

from concurrent.futures import ThreadPoolExecutor
import inspect
import math
from pathlib import Path

import libephemeris as eph
import pytest
from data_paths import DATA_ROOT


SATELLITE = DATA_ROOT / "tmp/planet_centers/plu017.bsp"
MAIN = DATA_ROOT / "de440.bsp"


@pytest.fixture(autouse=True)
def configured_engine():
    eph.close_moon_kernels()
    eph.set_calc_mode("skyfield")
    eph.set_precision_tier("medium")
    eph.set_leb_file(None)
    eph.set_ephe_path(str(MAIN.parent))
    eph.set_ephemeris_file(MAIN.name)
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_sid_mode(0)
    eph.set_topo(12.5, 41.9, 150.0)
    yield
    eph.close_moon_kernels()


def require_kernel():
    if not SATELLITE.is_file() or not MAIN.is_file():
        pytest.skip("local Pluto satellite and main DE440 kernels are required")


def test_ids_signatures_and_absent_target_do_not_open_a_main_image():
    assert eph.get_moon_name(9051) == eph.get_moon_name(9901) == "Charon"
    assert eph.get_moon_name(9501) == "Io"
    assert eph.get_moon_name(-1) == "Unknown Moon (-1)"
    assert eph.is_planetary_moon(9501) and not eph.is_planetary_moon(0)
    assert inspect.signature(eph.register_moon_spk).parameters["moons"].default is None
    assert list(inspect.signature(eph.calc_moon_position).parameters) == ["t", "moon_id", "iflag"]
    assert eph.list_registered_moons() == {}
    assert eph.get_moon_coverage(9901) is None
    eph.set_calc_mode("leb")
    assert eph.calc_moon_position(2451545.0, 9901, 258) is None


def test_registration_aliases_coverage_and_retirement():
    require_kernel()
    eph.register_moon_spk(str(SATELLITE))
    mapping = eph.list_registered_moons()
    assert {9051, 9901, 9902, 9903} <= mapping.keys()
    assert all(path == str(SATELLITE) for path in mapping.values())
    mapping.clear()
    assert eph.list_registered_moons()
    span = eph.get_moon_coverage(9901)
    assert span == eph.get_moon_coverage(9051)
    assert span[0] < 2451545.0 < span[1]
    eph.unregister_moon_spk(SATELLITE.name)
    assert eph.list_registered_moons() == {}
    eph.register_moon_spk(str(SATELLITE), [9051])
    assert set(eph.list_registered_moons()) == {9051, 9901}
    eph.close_moon_kernels()
    assert eph.get_moon_coverage(9901) is None


def test_registration_refuses_missing_files_and_images_without_selected_moons():
    with pytest.raises(FileNotFoundError):
        eph.register_moon_spk("not-an-available-satellite-file.bsp")
    require_kernel()
    with pytest.raises(ValueError):
        eph.register_moon_spk(str(SATELLITE), [9501])
    assert eph.list_registered_moons() == {}


def test_direct_helper_and_public_positions_use_registered_satellite():
    require_kernel()
    eph.register_moon_spk(str(SATELLITE))
    direct = eph.calc_moon_position(2451545.0, 9901, 258)
    assert isinstance(direct, tuple) and len(direct) == 6
    assert direct == eph.calc_moon_position(2451545.0, 9051, 258)
    assert all(math.isfinite(value) for value in direct)
    assert 0.0 <= direct[0] < 360.0 and -90.0 <= direct[1] <= 90.0
    assert direct[2] > 10.0
    for function in (eph.calc, eph.calc_ut):
        values, flags = function(2451545.0, 9901, 258)
        assert isinstance(values, tuple) and len(values) == 6
        assert all(math.isfinite(value) for value in values)
        assert isinstance(flags, int)
        pluto = function(2451545.0, 9, 258)[0]
        assert abs(values[0] - pluto[0]) < 0.01
        assert values != pluto


def test_barycentric_true_state_is_independent_of_the_observer_site():
    require_kernel()
    eph.register_moon_spk(str(SATELLITE))
    flags = 258 | 16 | 16384
    before = eph.calc_moon_position(2451545.0, 9901, flags)
    eph.set_topo(-120.0, -35.0, 2000.0)
    assert eph.calc_moon_position(2451545.0, 9901, flags) == before


def test_registered_moon_positions_are_deterministic_across_eight_workers():
    require_kernel()
    eph.register_moon_spk(str(SATELLITE))
    requests = [(2451545.0 + index, moon, 258) for index in range(4) for moon in (9901, 9902, 9903)]
    expected = [eph.calc_moon_position(*request) for request in requests]
    with ThreadPoolExecutor(max_workers=8) as pool:
        observed = list(pool.map(lambda request: eph.calc_moon_position(*request), requests))
    assert observed == expected
