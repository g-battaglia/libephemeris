"""Civil calendars, atomic offsets, leap seconds and native solar clocks."""
from data_paths import CORE_LEB

from pathlib import Path
import math

import pytest
import libephemeris as eph
from libephemeris import time_utils as time


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


@pytest.mark.parametrize("jd", [float("nan"), float("inf"), -float("inf")])
def test_calendar_rejects_nonfinite_epochs_before_integer_conversion(jd):
    with pytest.raises(ValueError):
        time.revjul(jd)


@pytest.mark.parametrize("jd", [1e300, -1e300])
def test_calendar_refuses_an_unrepresentable_civil_year(jd):
    with pytest.raises(OverflowError):
        time.revjul(jd)


def test_calendars_and_weekday():
    assert time.julday(2000, 1, 1) == 2451545.0
    assert time.revjul(2451545.0) == (2000, 1, 1, 12.0)
    assert time.julday(1582, 10, 5, cal=0) == time.julday(1582, 10, 15, cal=1)
    assert time.day_of_week(time.julday(2000, 1, 1)) == 5
    with pytest.raises(ValueError):
        time.julday(2000, 1, 1, cal=2)


def test_package_date_conversion_keeps_its_distinct_checked_layout():
    assert eph.date_conversion(2000, 1, 1) == (True, 2451545.0, (2000, 1, 1, 12.0))
    assert eph.date_conversion(2024, 2, 30, cal=b"g") == (False, time.julday(2024, 3, 1), (2024, 3, 1, 12.0))
    assert eph.date_conversion(1582, 10, 5, cal="J")[1] == eph.date_conversion(1582, 10, 15, cal="G")[1]
    with pytest.raises(ValueError):
        eph.date_conversion(2000, 1, 1, cal="gg")


def test_exact_atomic_offset_and_tai_helpers(configured_engine):
    assert time.TT_TAI_OFFSET_SECONDS == 32.184
    tt = time.tai_to_tt_jd(2451545.0)
    assert tt == 2451545.0 + 32.184 / 86400.0
    assert time.tt_to_tai_jd(tt) == 2451545.0
    utc = time.julday(2020, 1, 1, 0.0)
    assert time.get_tai_utc_for_jd(utc) == 37.0
    tai = time.utc_jd_to_tai(utc)
    assert time.tai_to_utc_jd(tai) == utc
    assert time.utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0) == tai
    assert time.tai_jd_to_utc(tai)[:5] == (2020, 1, 1, 0, 0)
    assert time.tai_jd_to_utc(tai)[5] == pytest.approx(0.0, abs=5e-5)


def test_utc_leap_second_and_independent_modern_clock(configured_engine):
    tt, ut1 = time.utc_to_jd(2016, 12, 31, 23, 59, 60.25)
    expected = time.jdet_to_utc(tt)
    assert expected[:5] == (2016, 12, 31, 23, 59)
    assert expected[5] == pytest.approx(60.25, abs=5e-5)
    assert time.jdut1_to_utc(ut1)[5] == pytest.approx(60.25, abs=5e-5)
    eph.set_delta_t_userdef(1000.0 / 86400.0)
    assert time.utc_to_jd(2016, 12, 31, 23, 59, 60.25) == (tt, ut1)
    with pytest.raises(eph.Error):
        time.utc_to_jd(2024, 1, 1, 23, 59, 60.0)


def test_timezone_changes_calendar_without_rounding_seconds():
    assert time.utc_time_zone(2024, 1, 1, 0, 15, 12.25, 1.0) == (2023, 12, 31, 23, 15, 12.25)
    assert time.utc_time_zone(2024, 1, 1, 12, 0, 0.0, -5.5) == (2024, 1, 1, 17, 30, 0.0)


def test_configured_delta_and_local_solar_clocks(configured_engine):
    delta = 64.0 / 86400.0
    assert time.deltat(2451545.0) == delta
    assert time.deltat_ex(2451545.0, flag=4) == delta
    gst = time.sidtime(2451545.0)
    assert 0.0 <= gst < 24.0
    assert time.sidtime(2451545.0, longitude=15.0) == pytest.approx((gst + 1.0) % 24.0, abs=1e-12)
    assert 0.0 <= time.sidtime0(2451545.0, 23.4393, 0.001) < 24.0
    equation = time.time_equ(2451545.0)
    assert abs(equation) < 0.02
    lmt = time.lat_to_lmt(2451545.0, 12.5)
    assert time.lmt_to_lat(lmt, 12.5) == pytest.approx(2451545.0, abs=1e-9)


def test_explicit_iau2006_reexports_native_mean_and_apparent_observables():
    assert eph.sidtime_iau2006 is time.sidtime_iau2006
    assert time.sidtime_iau2006(2451545., 2451545.) == pytest.approx(
        18.697138155539015, abs=1e-11)
    assert time.sidtime_iau2006(2451545., 2451545., apparent=False) == pytest.approx(
        18.69737482696563, abs=1e-11)
    with pytest.raises(TypeError):
        time.sidtime_iau2006(2451545., 2451545., False)


@pytest.mark.parametrize('apparent', [False, True])
def test_explicit_iau2006_clocks_are_independent_of_configured_delta_t(configured_engine, apparent):
    a = time.sidtime_iau2006(2451545., 2451545., apparent=apparent)
    eph.set_delta_t_userdef(1000. / 86400.)
    assert time.sidtime_iau2006(2451545., 2451545., apparent=apparent) == a
    b = time.sidtime_iau2006(2451546., 2451545., apparent=apparent)
    assert (b-a) % 24. == pytest.approx(24.*.00273781191135448, abs=1e-12)
    assert time.sidtime_iau2006(2451545., 2451546., apparent=apparent) != a


@pytest.mark.parametrize('ut1,tt', [
    (float('nan'), 2451545.), (2451545., float('inf')),
    (2451545., math.nextafter(2396757.5, -math.inf)),
    (2451545., math.nextafter(2469808.5, math.inf)),
    (math.nextafter(2451545.-73050000., -math.inf), 2451545.),
    (math.nextafter(2451545.+73050000., math.inf), 2451545.),
])
def test_explicit_iau2006_rejects_only_its_declared_checked_domain(ut1, tt):
    with pytest.raises(ValueError):
        time.sidtime_iau2006(ut1, tt)


@pytest.mark.parametrize('centuries', [-2001, 2001])
def test_supplied_angle_sidereal_preserves_finite_domain_beyond_origin_prototype(centuries):
    """The private origin's ±2000 guard is not a public API-domain contract."""
    eph.set_calc_mode('skyfield')
    eph.set_delta_t_model('espenak_meeus')
    eph.set_delta_t_userdef(None)
    eph.set_iers_delta_t_enabled(False)
    try:
        jd = 2451545.+centuries*36525.
        value = time.sidtime0(jd, 23.4, 0.)
        assert math.isfinite(value) and 0. <= value < 24.
    finally:
        eph.set_calc_mode(None)
        eph.set_delta_t_model(None)
        eph.set_delta_t_userdef(None)
        eph.set_iers_delta_t_enabled(None)


def test_legacy_finite_clamped_pole_domain_is_not_silently_replaced_by_an_origin_guard():
    """Computational admission is broader than unit-pole scientific validity."""
    eph.set_calc_mode('skyfield')
    eph.set_delta_t_userdef(0.)
    eph.set_iers_delta_t_enabled(False)
    try:
        jd = 2451545.+10000*36525.
        value = time.sidtime0(jd, 23.4, 0.)
        assert math.isfinite(value) and 0. <= value < 24.
    finally:
        eph.set_calc_mode(None)
        eph.set_delta_t_userdef(None)
        eph.set_iers_delta_t_enabled(None)
