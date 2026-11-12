"""Local minor-body identities, native resonance and heliocentric propagation."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path

import libephemeris as eph
from libephemeris import minor_bodies as bodies
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.clear_asteroid_name_cache()
    eph.clear_asteroid_elements_cache()


def test_local_names_and_cache_are_owned_by_the_native_engine():
    for name, number in ((" Ceres ", 1), ("Eris", 136199), ("chiron", 2060), ("Europa", 52)):
        assert eph.get_asteroid_number(name) == number
    assert eph.get_asteroid_number("") is None
    assert eph.clear_asteroid_name_cache() > 0
    assert eph.clear_asteroid_name_cache() == 0
    assert eph.clear_asteroid_elements_cache() == 0


def test_resonance_uses_the_shared_typed_minor_body_elements():
    elements = bodies.OrbitalElements("plutino", 2451545.0, 39.45, 0.2, 10.0, 40.0, 65.0, 30.0, 0.004)
    result = eph.detect_mean_motion_resonance(elements)
    assert result.is_resonant and result.resonance is not None
    assert (result.resonance.p, result.resonance.q) == (2, 3)
    assert result.resonance.a_resonant > 39.0
    assert elements.validate() == []


def test_major_catalog_and_non_downloadable_requests():
    rows = eph.list_major_asteroids()
    assert isinstance(rows, list) and rows
    for body, name in rows:
        assert eph.is_major_asteroid(body)
        info = eph.get_major_asteroid_info(body)
        assert info[3] == name
    assert eph.get_major_asteroid_info(-1) is None
    assert eph.auto_download_asteroid_spk(-1) is None
    assert eph.ensure_major_asteroid_spk(-1, 2451545.0) is False


def test_curated_numbered_asteroids_are_finite_and_thread_safe():
    requests = [(number, 2451545.0, False) for number in (1, 2, 3, 4, 2060, 5145)]
    expected = [eph.calc_asteroid_by_number(*request) for request in requests]
    for lon, lat, radius in expected:
        assert all(math.isfinite(value) for value in (lon, lat, radius))
        assert 0.0 <= lon < 360.0 and -90.0 <= lat <= 90.0 and radius > 0.0
    with ThreadPoolExecutor(max_workers=8) as pool:
        observed = list(pool.map(lambda request: eph.calc_asteroid_by_number(*request), requests))
    assert observed == expected
    with pytest.raises(ValueError):
        eph.calc_asteroid_by_number(0, 2451545.0)
