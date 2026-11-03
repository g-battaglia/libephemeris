"""Physical utility identities and the public numeric/text boundary."""
from data_paths import CORE_LEB

from pathlib import Path
import math

import pytest
import libephemeris as eph
from libephemeris import utils


@pytest.fixture
def configured_engine():
    path = CORE_LEB
    eph.set_precision_tier("base")
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(path))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    yield
    eph.set_leb_file(None)
    eph.set_delta_t_userdef(None)
    eph.set_iers_delta_t_enabled(None)
    eph.set_calc_mode(None)


def test_angle_wrap_and_midpoints():
    assert utils.degnorm(-90.0) == 270.0
    assert utils.radnorm(-math.pi) == math.pi
    assert utils.difdeg2n(350.0, 10.0) == -20.0
    assert utils.difdegn(350.0, 10.0) == 340.0
    assert utils.difrad2n(1.5 * math.pi, 0.5 * math.pi) == -math.pi
    assert utils.deg_midp(350.0, 10.0) == 0.0
    assert utils.rad_midp(1.5 * math.pi, 0.5 * math.pi) == 0.0


def test_integer_angles_and_display():
    assert utils.csnorm(-1) == 129599999
    assert utils.difcs2n(64800000, 0) == -64800000
    assert utils.difcsn(-1, 0) == 129599999
    assert utils.csroundsec(150) == 200
    assert utils.cs2degstr(3723456) == ' 10°20\'34.56"'
    assert utils.cs2lonlatstr(-360000, b"E", b"W") == "1W00"
    assert utils.cs2timestr(360000) == "01:00:00"
    assert utils.cs2timestr(360000, b"-", True) == "01-00"
    with pytest.raises(UnicodeDecodeError):
        utils.cs2timestr(0, b"\xff")
    assert utils.d2l(-1.5) == -2
    assert utils.split_deg(45.5) == (45, 30, 0, 0.0, 1)


def test_rotation_round_trip_and_rate_layouts():
    point = (123.4, -12.5, 2.0)
    rotated = utils.cotrans(point, 23.4393)
    restored = utils.cotrans(rotated, -23.4393)
    assert restored == pytest.approx(point, abs=1e-11)
    speed = (0.5, -0.1, 0.01)
    flat = utils.cotrans_sp((*point, *speed), 23.4393)
    split = utils.cotrans_sp(point, speed, 23.4393)
    assert flat == (*split[0], *split[1])
    assert flat[:3] == pytest.approx(rotated, abs=1e-11)
    step = 1e-5
    future = utils.cotrans(tuple(value + step * rate for value, rate in zip(point, speed)), 23.4393)
    for index in range(2):
        assert (future[index] - flat[index]) / step == pytest.approx(flat[index + 3], abs=1e-5)
    assert flat[2] == point[2] and flat[5] == speed[2]


def test_refraction_has_physical_direction_and_detail_units():
    assert 45.0 < utils.refrac(45.0) < 46.0
    assert 44.0 < utils.refrac(45.0, flag=1) < 45.0
    converted, details = utils.refrac_extended(20.0, 100.0)
    true, apparent, correction, dip = details
    assert converted == apparent
    assert true == 20.0
    assert apparent - true == pytest.approx(correction, abs=1e-14)
    assert correction > 0.0 and dip < 0.0


def test_horizontal_round_trip_and_angle_dictionary(configured_engine):
    geopos = (12.5, 41.9, 0.0)
    azimuth, altitude, _ = utils.azalt(2451545.0, 1, geopos, 1013.25, 15.0, (120.0, 10.0, 1.0))
    assert utils.azalt_rev(2451545.0, 1, geopos, azimuth, altitude) == pytest.approx((120.0, 10.0), abs=1e-10)
    angles = utils.calc_angles(2451545.0, geopos[1], geopos[0])
    assert len(angles) == 12
    assert angles["Descendant"] == pytest.approx((angles["Ascendant"] + 180.0) % 360.0)
    assert all(isinstance(value, float) and math.isfinite(value) for value in angles.values())


def test_checked_numeric_domains():
    with pytest.raises(ValueError):
        utils.d2l(math.nan)
    with pytest.raises(OverflowError):
        utils.d2l(math.inf)
    with pytest.raises(ValueError):
        utils.split_deg(math.nan)
    with pytest.raises(OverflowError):
        utils.csnorm(2**100)
