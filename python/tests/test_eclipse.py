"""Native event records, physical contact order and public tuple layouts."""
from data_paths import CORE_LEB

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


def test_besselian_record_and_native_interpolation():
    record = eph.BesselianElements(2451545.0, 0.1, 0.2, 5.0, 0.6, -0.01, 359.0,
        0.1, 0.2, 0.03, 0.04, 0.005, 2.0)
    advanced = eph.interpolate_besselian_elements(record, 2451545.125)
    assert advanced.t0 == 2451545.125
    assert advanced.x == pytest.approx(0.4)
    assert advanced.y == pytest.approx(0.8)
    assert advanced.mu == 5.0
    assert advanced.dx_dt == record.dx_dt
    record.x = 0.5
    assert record.x == 0.5
    assert list(inspect.signature(eph.BesselianElements).parameters) == [
        "t0", "x", "y", "d", "l1", "l2", "mu", "dx_dt", "dy_dt", "dd_dt", "dl1_dt", "dl2_dt", "dmu_dt"]


def test_fixed_instant_event_layouts():
    flag, geography, attributes = eph.sol_eclipse_where(2460409.264)
    assert isinstance(flag, int)
    assert isinstance(geography, tuple) and len(geography) == 10
    assert isinstance(attributes, tuple) and len(attributes) == 20
    assert all(math.isfinite(value) for value in geography + attributes)
    assert geography[2:] == (0.0,) * 8
    lunar_flag, lunar = eph.lun_eclipse_how(2459891.958, (0.0, 0.0, 0.0))
    assert isinstance(lunar_flag, int) and isinstance(lunar, tuple) and len(lunar) == 20
    assert lunar[0] == lunar[8]


def test_lunar_contact_order_and_duration_units():
    maximum = 2459891.958
    p1 = eph.calc_lunar_eclipse_penumbral_first_contact_p1(maximum)
    u1 = eph.calc_lunar_eclipse_umbral_first_contact_u1(maximum)
    u2 = eph.calc_lunar_eclipse_umbral_second_contact_u2(maximum)
    u3 = eph.calc_lunar_eclipse_umbral_third_contact_u3(maximum)
    u4 = eph.calc_lunar_eclipse_umbral_fourth_contact_u4(maximum)
    p4 = eph.calc_lunar_eclipse_penumbral_fourth_contact_p4(maximum)
    assert p1 < u1 < u2 < maximum < u3 < u4 < p4
    assert eph.calc_lunar_eclipse_total_duration(maximum) == pytest.approx((u3 - u2) * 24 * 60, abs=1e-8)
    assert eph.calc_lunar_eclipse_umbral_duration(maximum) == pytest.approx((u4 - u1) * 24 * 60, abs=1e-8)


def test_rise_set_and_named_star_share_public_layout():
    for body in (0, "Sirius"):
        flag, times = eph.rise_trans(2460389.5, body, 4, (12.5, 41.9, 0.0))
        assert flag == 0 and isinstance(times, tuple) and len(times) == 10
        assert 2460389.5 <= times[0] < 2460390.5
        assert times[1:] == (0.0,) * 9


def test_event_signatures_and_invalid_domains():
    assert inspect.signature(eph.sol_eclipse_when_glob).parameters["backwards"].default is False
    assert inspect.signature(eph.lun_occult_where).parameters["body"].default == 0
    with pytest.raises(ValueError):
        eph.rise_trans(2451545.0, 0, 1, (12.5, 41.9))
    with pytest.raises(ValueError):
        eph.get_saros_number(2451545.0, "invalid")
    with pytest.raises(ValueError):
        eph.sol_eclipse_max_time(2460409.264, lat=32.8)
    with pytest.raises(ValueError):
        eph.planet_occult_when_glob(2451545.0, 3)
    with pytest.raises(eph.Error):
        eph.lun_occult_where(2451545.0, 1)
