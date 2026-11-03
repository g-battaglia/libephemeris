"""Independent contexts keep settings and snapshots through concurrent work."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import inspect

import libephemeris as eph
import pytest

FIXTURE = CORE_LEB


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(FIXTURE))
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_iers_delta_t_enabled(False)


def test_independent_settings_and_cache():
    left, right = eph.EphemerisContext(), eph.EphemerisContext()
    left.set_topo(12.0, 42.0, 100.0)
    assert left.get_topo() == (12.0, 42.0, 100.0)
    assert right.get_topo() is None
    left.set_sid_mode(255, 2451545.0, 24.0)
    assert left.get_sid_mode(True) == (255, 2451545.0, 24.0)
    assert right.get_sid_mode() == 0
    values = {"ascendant": 12.5}
    left.set_angles_cache(values)
    values["ascendant"] = 99.0
    assert left.get_angles_cache() == {"ascendant": 12.5}
    assert right.get_angles_cache() == {}
    left.clear_angles_cache()
    assert left.get_angles_cache() == {}


def test_context_calls_equal_global_and_survive_class_close():
    context = eph.EphemerisContext()
    requests = [(2451545.0 + day, body, 258) for day in range(3) for body in range(10)]
    expected = [eph.calc_ut(*row) for row in requests]
    assert [context.calc_ut(*row) for row in requests] == expected
    with ThreadPoolExecutor(max_workers=8) as workers:
        actual = list(workers.map(lambda row: context.calc_ut(*row), requests))
    assert actual == expected
    context.set_topo(12.0, 42.0)
    retained = context.get_leb_reader()
    eph.EphemerisContext.close()
    assert context.get_topo() == (12.0, 42.0, 0.0)
    assert retained.has_body(0)
    assert context.calc_ut(*requests[0]) == expected[0]


def test_native_handle_shapes_and_method_signatures():
    context = eph.EphemerisContext()
    assert isinstance(context.get_loader().directory, str)
    whole, fraction = context.get_timescale().ut1_jd(2451545.0)
    assert whole == 2451545.0 and 0 < fraction < 0.01
    assert context.get_timescale().tt_jd(2451545.125) == (2451545.0, 0.125)
    assert list(inspect.signature(eph.EphemerisContext).parameters) == ["ephe_path", "ephe_file"]
    assert list(inspect.signature(context.calc_ut).parameters) == ["tjd_ut", "ipl", "iflag"]
    assert inspect.signature(context.get_sid_mode).parameters["full"].default is False
    assert context.get_spk_body_info(99999) is None
    context.unregister_spk_body(99999)
    assert context.list_spk_bodies() == {}


def test_invalid_coordinates_do_not_replace_observer():
    context = eph.EphemerisContext()
    context.set_topo(12.0, 42.0)
    with pytest.raises(eph.CoordinateError):
        context.set_topo(12.0, 91.0)
    assert context.get_topo() == (12.0, 42.0, 0.0)
