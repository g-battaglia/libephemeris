"""Geometric house identities, polar fallback and public result layouts."""
from data_paths import CORE_LEB

import importlib
import inspect
from pathlib import Path

import libephemeris as eph
import pytest

house = importlib.import_module("libephemeris.houses")


@pytest.fixture(autouse=True)
def configured_engine():
    eph.set_precision_tier("base")
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)
    eph.set_sid_mode(0)


def test_root_function_and_module_coexist_with_original_defaults():
    assert callable(eph.houses)
    assert eph.houses is house.houses
    assert inspect.signature(house.houses).parameters["hsys"].default == ord("P")
    assert inspect.signature(house.houses_with_fallback).parameters["fallback_hsys"].default == ord("O")
    assert inspect.signature(house.gauquelin_sector).parameters["flags"].default == 32770


def test_equal_and_whole_sign_house_geometry_and_layouts():
    cusps, angles = house.houses_armc(123.45, 41.9, 23.44, b"E")
    assert isinstance(cusps, tuple) and isinstance(angles, tuple)
    assert (len(cusps), len(angles)) == (12, 8)
    assert cusps[0] == angles[0]
    for first, second in zip(cusps, (*cusps[1:], cusps[0])):
        assert (second - first) % 360.0 == pytest.approx(30.0, abs=1e-12)
    position = house.house_pos(123.45, 41.9, 23.44, (cusps[0] + 15.0, 0.0), "E")
    assert position == pytest.approx(1.5, abs=1e-12)
    assert position == house.house_pos(123.45, 41.9, 23.44, ord("E"), cusps[0] + 15.0)
    whole, _ = house.houses_armc(123.45, 41.9, 23.44, "W")
    assert all(value % 30.0 == 0.0 for value in whole)
    sectors, _ = house.houses_armc(123.45, 41.9, 23.44, "G")
    assert len(sectors) == 36


def test_extended_positions_and_rate_layouts():
    plain = house.houses(2451545.0, 41.9, 12.5)
    extended = house.houses_ex(2451545.0, 41.9, 12.5)
    rates = house.houses_ex2(2451545.0, 41.9, 12.5)
    assert plain == extended == rates[:2]
    assert tuple(map(len, rates)) == (12, 8, 12, 8)
    assert all(isinstance(part, tuple) for part in rates)
    armc_rates = house.houses_armc_ex2(123.45, 41.9, 23.44)
    assert armc_rates[:2] == house.houses_armc(123.45, 41.9, 23.44)
    assert tuple(map(len, armc_rates)) == (12, 8, 12, 8)


def test_polar_fallback_and_latitude_information():
    with pytest.raises(eph.PolarCircleError):
        house.houses_armc(123.45, 70.0, 23.44, "P")
    fallback = house.houses_armc_with_fallback(123.45, 70.0, 23.44, "P")
    assert fallback[:2] == house.houses_armc(123.45, 70.0, 23.44, "O")
    assert fallback[2] is True and "Porphyry" in fallback[3]
    fallback = house.houses_with_fallback(2451545.0, 70.0, 12.5, "P")
    assert fallback[:2] == house.houses(2451545.0, 70.0, 12.5, "O")
    assert fallback[2] is True
    assert house.get_polar_latitude_threshold() == pytest.approx(66.56)
    info = house.get_extreme_latitude_info(-85.0)
    assert info["hemisphere"] == "S"
    assert info["is_extreme"] and info["is_polar_circle"]
    assert info["affected_systems"] == ["P", "K", "G"]
    assert house.house_name(b"P") == "Placidus"


@pytest.mark.parametrize("api", [house.houses_armc, house.houses_armc_ex2])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
@pytest.mark.parametrize("latitude,epsilon", [
    (float.fromhex("0x1.0a3e2a79504e0p+6"), 23.4392911), (66.56, 23.44),
])
def test_polar_metadata_and_refusal_preserve_the_exact_input_side(api, sign, latitude, epsilon):
    assert latitude + epsilon == 90.0  # Rounded addition loses the strict exterior.
    info = house.get_extreme_latitude_info(sign * latitude, epsilon)
    assert info["is_polar_circle"] and info["affected_systems"] == ["P", "K", "G"]
    with pytest.raises(eph.PolarCircleError):
        api(90.0, sign * latitude, epsilon, "P")


@pytest.mark.parametrize("api", [house.houses_armc, house.houses_armc_ex2])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
@pytest.mark.parametrize("latitude,epsilon", [
    (float.fromhex("0x1.0a3e2a79504dfp+6"), 23.4392911), (66.5, 23.5),
])
def test_equatorward_neighbour_and_exact_equality_remain_admissible(api, sign, latitude, epsilon):
    assert not house.get_extreme_latitude_info(sign * latitude, epsilon)["is_polar_circle"]
    assert len(api(90.0, sign * latitude, epsilon, "P")[0]) == 12


def test_native_gauquelin_sectors_are_within_a_diurnal_turn():
    for body in (0, "Spica"):
        for method in (0, 1):
            sector = house.gauquelin_sector(2451545.0, body, method, (12.5, 41.9, 0.0))
            assert 1.0 <= sector < 37.0


def test_sunshine_polar_fallback_is_the_exact_porphyry_definition():
    for latitude in (-85.0, -70.0, 70.0, 85.0):
        fallback = house.houses_with_fallback(2451545.0, latitude, 12.5, "i")
        assert fallback[:2] == house.houses(2451545.0, latitude, 12.5, "O")
        assert fallback[2] is True and "Porphyry" in fallback[3]
        cusps, angles = fallback[:2]
        assert cusps[0] == angles[0] and cusps[9] == angles[1]
        for quadrant in range(4):
            first = cusps[3 * quadrant]
            last = cusps[(3 * quadrant + 3) % 12]
            third = (last - first) % 360.0 / 3.0
            for division in (1, 2):
                assert cusps[3 * quadrant + division] == pytest.approx((first + third * division) % 360.0, abs=1e-12)


def test_dynamic_sidereal_houses_share_the_native_ayanamsha():
    jd, latitude, longitude = 2460310.5, -33.9, 151.2
    tropical, angles = house.houses_ex(jd, latitude, longitude, "P")
    for mode in (27, 30):
        eph.set_sid_mode(mode)
        offset = eph.get_ayanamsa_ex_ut(jd, 0)[1]
        sidereal, shifted_angles = house.houses_ex(jd, latitude, longitude, "P", 65536)
        assert shifted_angles[2] == angles[2]
        for before, after in zip(tropical, sidereal):
            assert after == pytest.approx((before - offset) % 360.0, abs=1e-12)
