"""Lunar geometric and analytic invariants at scalar TT epochs."""
from data_paths import CORE_LEB

import math
from pathlib import Path

import libephemeris as eph
from libephemeris import lunar
import pytest


@pytest.fixture(autouse=True)
def configured_engine():
    eph.set_precision_tier("base")
    eph.set_calc_mode("leb")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)


def test_lunar_warning_and_error_types():
    assert issubclass(lunar.MeeusPolynomialWarning, UserWarning)
    assert issubclass(lunar.MeeusRangeError, ValueError)


def test_mean_node_regresses_and_mean_apogee_is_an_inclined_point():
    first = lunar.calc_mean_lunar_node(2451544.5)
    second = lunar.calc_mean_lunar_node(2451545.5)
    rate = (second - first + 180.0) % 360.0 - 180.0
    assert -0.06 < rate < -0.04
    longitude, latitude = lunar.calc_mean_lilith_with_latitude(2451545.0)
    assert longitude == lunar.calc_mean_lilith(2451545.0)
    assert 0.0 <= longitude < 360.0 and abs(latitude) < 6.0
    with pytest.raises(eph.EphemerisRangeError):
        lunar.calc_mean_lunar_node(math.nan)


def test_true_node_is_on_ecliptic_and_apogee_is_the_conic_maximum():
    node = lunar.calc_true_lunar_node(2451545.0)
    apogee = lunar.calc_true_lilith(2451545.0)
    assert node[1] == 0.0
    assert 0.002 < node[2] < 0.003
    assert node[2] <= apogee[2]
    assert 0.0 <= node[0] < 360.0 and 0.0 <= apogee[0] < 360.0
    assert len(node) == len(apogee) == 3


def test_interpolated_apses_have_distinct_physical_distances():
    apogee = lunar.calc_interpolated_apogee(2451545.0)
    perigee = lunar.calc_interpolated_perigee(2451545.0)
    assert apogee[2] > perigee[2] > 0.0
    assert all(math.isfinite(value) for value in (*apogee, *perigee))
    assert 0.0 <= apogee[0] < 360.0 and 0.0 <= perigee[0] < 360.0


def test_tt_arguments_are_independent_of_configured_ut1_offset():
    functions = (lunar.calc_true_lunar_node, lunar.calc_true_lilith,
                 lunar.calc_mean_lunar_node, lunar.calc_mean_lilith,
                 lunar.calc_mean_lilith_with_latitude,
                 lunar.calc_interpolated_apogee, lunar.calc_interpolated_perigee)
    expected = [function(2451545.0) for function in functions]
    eph.set_delta_t_userdef(0.02)
    assert [function(2451545.0) for function in functions] == expected
