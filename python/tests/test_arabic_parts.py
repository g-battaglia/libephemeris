"""Classical lots, missing inputs and geometric sect at polar latitudes."""
from data_paths import CORE_LEB

from pathlib import Path
import libephemeris as eph
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)


def positions(sun):
    return {"Asc": 0.0, "Sun": sun, "Moon": 240.0, "Mercury": 130.0, "Venus": 110.0}


def test_lots_reverse_fortune_and_spirit_with_sect():
    day = eph.calc_all_arabic_parts(positions(270.0))
    assert day == {"Pars_Fortunae": 330.0, "Pars_Spiritus": 30.0, "Pars_Amoris": 200.0, "Pars_Fidei": 250.0}
    night = eph.calc_all_arabic_parts(positions(90.0))
    assert night["Pars_Fortunae"] == 210.0 and night["Pars_Spiritus"] == 150.0
    for sun in (0.0, 180.0):
        result = eph.calc_all_arabic_parts(positions(sun))
        assert result["Pars_Fortunae"] == (240.0 - sun) % 360.0


def test_high_latitude_sect_uses_the_true_solar_altitude():
    for latitude in (-85.0, -70.0, 70.0, 85.0):
        for sun in (0.0, 90.0, 180.0, 270.0):
            sample = positions(sun)
            result = eph.calc_all_arabic_parts(sample, jd=2460310.5, geo_lat=latitude, geo_lon=12.5)
            altitude = eph.azalt(2460310.5, 0, [12.5, latitude, 0.0], 0.0, 0.0, [sun, 0.0, 1.0])[1]
            expected = (240.0 - sun) % 360.0 if altitude >= 0.0 else (sun - 240.0) % 360.0
            assert result["Pars_Fortunae"] == expected


def test_missing_longitudes_are_errors():
    with pytest.raises(KeyError):
        eph.calc_all_arabic_parts({"Asc": 0.0, "Sun": 90.0})
