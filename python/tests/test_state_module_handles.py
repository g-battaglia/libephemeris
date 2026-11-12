"""State globals retain native resources and independent scalar configuration."""

from data_paths import CORE_LEB, DATA_ROOT

from concurrent.futures import ThreadPoolExecutor
import hashlib
import inspect
import json
from pathlib import Path
import subprocess

import libephemeris as e
import libephemeris.state as state
import pytest

FIXTURE = CORE_LEB


@pytest.fixture(autouse=True)
def configured(monkeypatch, tmp_path):
    monkeypatch.setenv("LIBEPHEMERIS_DATA_DIR", str(tmp_path))
    monkeypatch.setenv("LIBEPHEMERIS_CONFIG", str(tmp_path / "absent.toml"))
    monkeypatch.setenv("LIBEPHEMERIS_ENV_FILE", str(tmp_path / "absent.env"))
    e.close()
    e.reload_configuration()
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    e.set_leb_file(str(FIXTURE))
    yield
    monkeypatch.undo()
    e.close()
    e.reload_configuration()
    e.set_calc_mode("leb")
    e.set_precision_tier("base")
    e.set_leb_file(str(FIXTURE))


def test_state_scalar_getters_are_owned_and_contexts_remain_independent():
    independent = e.EphemerisContext()
    assert state.get_topo() is None
    e.set_topo(12.5, 41.9, 123.0)
    assert state.get_topo() == (12.5, 41.9, 123.0)
    assert independent.get_topo() is None
    e.set_sid_mode(e.SIDM_USER, 2451545.0, 12.0)
    assert state.get_sid_mode() == e.SIDM_USER
    assert state.get_sid_mode(True) == (e.SIDM_USER, 2451545.0, 12.0)
    values = {"Ascendant": 123.4, "Sun": 56.7}
    state.set_angles_cache(values)
    values["Sun"] = 0.0
    retained = state.get_angles_cache()
    assert retained == {"Ascendant": 123.4, "Sun": 56.7}
    retained["Sun"] = 99.0
    assert state.get_angles_cache()["Sun"] == 56.7
    independent.set_angles_cache({"Moon": 200.0})
    assert "Moon" not in state.get_angles_cache()
    with ThreadPoolExecutor(max_workers=8) as pool:
        assert (
            list(pool.map(lambda _: state.get_angles_cache(), range(24)))
            == [state.get_angles_cache()] * 24
        )
    state.clear_angles_cache()
    assert state.get_angles_cache() == {}
    assert independent.get_angles_cache() == {"Moon": 200.0}


def test_default_handles_use_real_configuration_and_sealed_main_policy(tmp_path):
    loader = state.get_loader()
    clock = state.get_timescale()
    assert loader.directory == str(tmp_path)
    assert clock.tt_jd(2451545.125) == (2451545.0, 0.125)
    whole, fraction = clock.ut1_jd(2451545.0)
    assert whole == 2451545.0 and 0.0 < fraction < 0.01
    assert state.get_planet_centers() is None
    assert state.get_planet_center_segment(599) is None
    assert not state.jpl_fallback_active()
    with pytest.raises(RuntimeError):
        state.get_planets()
    e.close()
    assert clock.tt_jd(2451545.125) == (2451545.0, 0.125)
    assert loader.directory == str(tmp_path)


def test_fallback_predicate_does_not_open_or_download_selected_files(tmp_path):
    e.set_calc_mode("skyfield")
    e.set_network_policy("sealed")
    e.set_ephemeris_file("absent-custom-kernel.bsp")
    assert state.jpl_fallback_active()
    assert list(tmp_path.iterdir()) == []
    e.set_ephemeris_file("de441.bsp")
    assert not state.jpl_fallback_active()
    explicit = tmp_path / "custom-kernel.bsp"
    explicit.write_bytes(b"format validation belongs to the opener")
    e.set_ephemeris_file(str(explicit))
    assert not state.jpl_fallback_active()
    assert explicit.read_bytes() == b"format validation belongs to the opener"
