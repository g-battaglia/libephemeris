"""Convenience helpers call native kernels and retain their public layouts."""
from data_paths import CORE_LEB

import math
import time
from pathlib import Path
import libephemeris as eph
from libephemeris import contrib
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)


def test_scalar_partitions_names_and_aspect_layouts():
    assert contrib.long2rasi(361.0) == 0
    assert contrib.long2nakshatra(0.0) == (0,0)
    assert contrib.rasi_dif(11,0) == 1
    assert contrib.rasi_dif2(0,11) == -1
    assert contrib.signtostr(0) == "Aries"
    assert contrib.house_system_name("SS") == "Unknown"
    assert contrib.match_aspect(0,0.0,90.0,1,90.5,1.0) == (1,0.5)
    assert contrib.match_aspect2(0.0,90.5,90.0,1.0) == (1,0.5)
    assert contrib.match_aspect3(0,0.0,90.0,1,90.5,1.0) == (1,0.5,90.5)
    assert contrib.match_aspect4(0,0.0,90.0,1,90.5,1.0) == (1,0.5,90.5,90.5)
    assert contrib.ASWINI == 0 and contrib.REVATHI == 26


def test_time_and_coordinate_formats():
    assert contrib.t2i("12") == (12,0,0)
    assert contrib.dt2i("2024-03-20T12:34:56.123456+02:00") == (2024,3,20,12,34,56)
    assert contrib.jd2isostr(2451545.0) == "2000-01-01T12:00:00"
    assert contrib.geoc2d("N30") == 0.5
    assert contrib.geoc2d("12S30") == -12.5
    assert contrib.geolat2c(-12.5) == "12S30"
    assert contrib.geolon2c(151.2) == "151E12"
    assert contrib.years_diff("2024-01-01T12:00:00+02:00", "2024-01-01T10:00:00Z") == 0.0
    assert contrib.years_diff("2024-01-01T00:00:00.100000", "2024-01-01T00:00:00.200000") * (365.25*86400) == pytest.approx(0.1,abs=1e-12)


def test_houses_stars_and_sampled_walks_use_shared_engine():
    assert contrib.raman_houses(2451545.0,12.5,41.9) == eph.houses(2451545.0,41.9,12.5,ord("S"))[0]
    stars = contrib.saturn_4_stars(2451545.0)
    assert isinstance(stars,tuple) and len(stars)==4 and all(math.isfinite(value) for value in stars)
    for name in ("next_aspect","next_aspect2","next_aspect_cusp","next_aspect_cusp2"):
        jd,residual=getattr(contrib,name)(2,90.0,0.0,2451545.0)
        assert jd>=2451545.0 and residual<=1.0
    jd,residual=contrib.next_aspect_with(2,3,0.0,2451545.0)
    assert jd>=2451545.0 and residual<=1.0
    assert contrib.next_aspect_with2(2,3,0.0,2451545.0)==(jd,residual)
    assert contrib.next_retro(2,2451545.0)>2451545.0


def test_clock_is_a_live_utc_julian_day_and_containers_are_plain():
    before=time.time()/86400+2440587.5
    result=contrib.jdnow()
    after=time.time()/86400+2440587.5
    assert before<=result<=after
    user=contrib.User();user.label="example"
    data=contrib.Data();data.value=1.0
    assert user.label=="example" and data.value==1.0
