"""The remaining planetary calls retain typed results and owned configuration."""

from data_paths import CORE_LEB

import inspect
import json
import math
from pathlib import Path

import libephemeris as ephem
from libephemeris import planets
import pytest


@pytest.fixture(autouse=True)
def explicit_context():
    ephem.set_calc_mode("leb")
    ephem.set_precision_tier("base")
    ephem.set_leb_file(str(CORE_LEB))
    ephem.set_iers_delta_t_enabled(False)
    ephem.set_delta_t_userdef(64.0 / 86400.0)
    ephem.set_sid_mode(0)
    yield
    ephem.set_sid_mode(0)


def test_tuple_layouts_and_native_values():
    values, flag = planets.calc_pctr(2451545.0, ephem.MOON, ephem.EARTH)
    assert type(values) is tuple and len(values) == 6 and type(flag) is int
    assert all(type(value) is float and math.isfinite(value) for value in values)
    for function in (planets.nod_aps, planets.nod_aps_ut):
        points = function(2451545.0, ephem.MARS)
        assert type(points) is tuple and len(points) == 4
        assert all(type(point) is tuple and len(point) == 6 for point in points)
        assert all(math.isfinite(value) for point in points for value in point)
    for function in (planets.get_orbital_elements, planets.get_orbital_elements_ut):
        elements = function(2451545.0, ephem.MARS, ephem.FLG_HELCTR)
        assert type(elements) is tuple and len(elements) == 50
        assert all(type(value) is float and math.isfinite(value) for value in elements)
        assert elements[17:] == (0.0,) * 33
    distances = planets.orbit_max_min_true_distance(
        2451545.0, ephem.MARS, ephem.FLG_HELCTR
    )
    assert type(distances) is tuple and len(distances) == 3
    assert distances[0] >= distances[1] > 0.0
    for function in (planets.pheno, planets.pheno_ut):
        values = function(2451545.0, ephem.MARS)
        assert type(values) is tuple and len(values) == 20
        assert 0.0 <= values[1] <= 1.0
        assert all(type(value) is float and math.isfinite(value) for value in values)


def test_elongation_helpers_share_the_native_calculation():
    value, evening = planets.get_elongation_from_sun(2451545.0, ephem.VENUS)
    assert type(value) is float and type(evening) is bool
    assert planets.get_signed_elongation(2451545.0, ephem.VENUS) == value
    assert planets.is_evening_star(2451545.0, ephem.VENUS) == evening
    assert planets.is_morning_star(2451545.0, ephem.VENUS) != evening
    assert planets.get_elongation_type(2451545.0, ephem.VENUS) == (
        "eastern" if evening else "western"
    )
    assert planets.is_morning_star(2451545.0, ephem.SUN) is False
    assert planets.is_evening_star(2451545.0, ephem.SUN) is False
    assert planets.get_elongation_type(2451545.0, ephem.SUN) == "none"


def test_sidereal_selection_changes_the_shared_runtime():
    planets.set_sid_mode(ephem.SIDM_USER, 2451545.0, 24.0)
    assert abs(planets.get_ayanamsa(2451545.0) - 24.0) < 1e-5
    assert type(planets.get_ayanamsa_ut(2451545.0)) is float
    for function in (planets.get_ayanamsa_ex, planets.get_ayanamsa_ex_ut):
        result = function(2451545.0)
        assert type(result) is tuple and len(result) == 2
        assert type(result[0]) is int and type(result[1]) is float
    planets.set_sid_mode(ephem.SIDM_LAHIRI)
    assert planets.get_ayanamsa_name(ephem.SIDM_LAHIRI) == "Lahiri"
    assert planets.get_ayanamsa_name(255) == ""


def test_native_names_and_backend_metadata():
    assert planets.get_planet_name(ephem.SUN) == "Sun"
    assert planets.get_planet_name(ephem.ERIS) == "Eris"
    assert planets.get_planet_name(ephem.MOON_IO) == "Io"
    assert planets.get_planet_name(ephem.PLMOON_IO) == "Io"
    assert planets.get_planet_name(-999) == ""
    assert planets.get_nutation_model() == {
        "model": "IAU2006_2000A",
        "terms": 1365,
        "precision": "~0.01-0.05 mas",
        "source": "native",
        "skyfield_available": False,
    }
    for warning in (planets.NutationFallbackWarning, planets.HeliocentricNodApsWarning):
        assert issubclass(warning, UserWarning)
        assert str(warning("notice")) == "notice"


def test_standard_argument_and_orbital_errors():
    with pytest.raises(TypeError):
        planets.pheno(object(), ephem.MARS)
    with pytest.raises(ephem.Error):
        planets.get_orbital_elements(2451545.0, ephem.SUN, 0)
    with pytest.raises(ephem.Error):
        planets.get_orbital_elements_ut(2451545.0, ephem.SUN, 0)
