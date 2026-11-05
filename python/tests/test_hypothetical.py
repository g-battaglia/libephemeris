"""Historical orbits, typed records, Keplerian invariants and aliases."""
from data_paths import CORE_LEB

from concurrent.futures import ThreadPoolExecutor
import inspect
import math
from pathlib import Path

import libephemeris as eph
from libephemeris import hypothetical as models
import pytest


@pytest.fixture(autouse=True)
def configured():
    eph.set_calc_mode("leb")
    eph.set_precision_tier("base")
    eph.set_leb_file(str(CORE_LEB))
    eph.set_iers_delta_t_enabled(False)
    eph.set_delta_t_userdef(64.0 / 86400.0)


def test_published_element_records_are_native_snapshots():
    assert set(models.URANIAN_KEPLERIAN_ELEMENTS) == set(range(40, 48))
    assert models.TRANSPLUTO_KEPLERIAN_ELEMENTS.a == 77.755
    assert models.VULCAN_ELEMENTS.Omega0 == 47.787931
    assert models.WALDEMATH_ELEMENTS.M_rate == 0.0
    assert isinstance(models.TRANSPLUTO_KEPLERIAN_ELEMENTS, models.TransplutoKeplerianElements)
    with pytest.raises(AttributeError):
        models.TRANSPLUTO_KEPLERIAN_ELEMENTS.a = 1.0


def test_polynomial_and_circular_orbit_have_exact_geometry():
    polynomial = models.TPolynomial(10.0, 2.0)
    assert polynomial.evaluate(5.0) == 20.0
    orbit = models.OrbitalElements("circular", 2451545.0, 2451545.0, False,
        models.TPolynomial(30.0), 2.0, models.TPolynomial(), models.TPolynomial(),
        models.TPolynomial(), models.TPolynomial())
    assert orbit.get_mean_motion() > 0.0
    lon, lat, radius, dlon, dlat, dr = models.calc_orbital_position(orbit, orbit.epoch_jd)
    assert lon == pytest.approx(30.0, abs=1e-12)
    assert lat == 0.0 and radius == pytest.approx(2.0)
    assert dlon == pytest.approx(orbit.get_mean_motion(), rel=1e-9)
    assert dlat == pytest.approx(0.0, abs=1e-12) and dr == pytest.approx(0.0, abs=1e-12)
    assert models.SeorbelElements is models.OrbitalElements


@pytest.mark.parametrize("axis,error", [
    (0.0, ZeroDivisionError),
    (-0.0, ZeroDivisionError),
    (float.fromhex("0x1p-720"), ZeroDivisionError),
    (float.fromhex("0x0.0000000000001p-1022"), ZeroDivisionError),
    (float.fromhex("0x1p+684"), OverflowError),
    (float.fromhex("0x1.fffffffffffffp+1023"), OverflowError),
], ids=["positive-zero", "negative-zero", "power-underflow", "subnormal-axis",
        "power-overflow", "maximum-axis"])
def test_mean_motion_reports_power_failures_on_construction_and_mutation(axis, error):
    orbit = models.OrbitalElements("power-domain", 2451545.0, 2451545.0, False,
        models.TPolynomial(30.0), axis, models.TPolynomial(), models.TPolynomial(),
        models.TPolynomial(), models.TPolynomial())
    assert orbit.semi_axis.hex() == axis.hex()
    with pytest.raises(error):
        orbit.get_mean_motion()
    orbit.semi_axis = 4.0
    ordinary = orbit.get_mean_motion()
    assert math.isfinite(ordinary) and ordinary > 0.0
    orbit.semi_axis = axis
    with pytest.raises(error):
        orbit.get_mean_motion()
    orbit.semi_axis = 4.0
    assert orbit.get_mean_motion().hex() == ordinary.hex()


@pytest.mark.parametrize("axis", [float.fromhex("0x1p-716"), float.fromhex("0x1p-700")])
def test_mean_motion_keeps_division_overflow_with_nonzero_power(axis):
    orbit = models.OrbitalElements("division-range", 2451545.0, 2451545.0, False,
        models.TPolynomial(), axis, models.TPolynomial(), models.TPolynomial(),
        models.TPolynomial(), models.TPolynomial())
    assert orbit.get_mean_motion() == math.inf


def test_mean_motion_preserves_ordinary_words_and_nonfinite_results():
    orbit = models.OrbitalElements("unchanged-domain", 2451545.0, 2451545.0, False,
        models.TPolynomial(), 1.0, models.TPolynomial(), models.TPolynomial(),
        models.TPolynomial(), models.TPolynomial())
    # Repeated evaluation and mutation must preserve this engine's own words;
    # libm power results are not specified by a platform-independent hex snapshot.
    for axis in (1.0, 2.0, 4.0, 16.0, 0.01, 77.755, 1e-200,
                 float.fromhex("0x1p+680")):
        orbit.semi_axis = axis
        expected = orbit.get_mean_motion()
        assert math.isfinite(expected) and expected > 0.0
        assert orbit.get_mean_motion().hex() == expected.hex()
        orbit.semi_axis = 4.0 if axis == 1.0 else 1.0
        assert math.isfinite(orbit.get_mean_motion())
        orbit.semi_axis = axis
        assert orbit.get_mean_motion().hex() == expected.hex()
    for axis in (math.inf, -math.inf):
        orbit.semi_axis = axis
        value = orbit.get_mean_motion()
        assert value == 0.0 and math.copysign(1.0, value) == 1.0
    for axis in (math.nan, -2.0):
        orbit.semi_axis = axis
        assert math.isnan(orbit.get_mean_motion())


def test_installed_bundled_csv_and_legacy_aliases():
    path = models.get_bundled_fictitious_orbits_path()
    assert isinstance(path, Path) and path.is_file() and path.name == "fictitious_orbits.csv"
    rows = models.load_bundled_fictitious_orbits()
    assert len(rows) == 12 and rows[0].name == "Cupido"
    assert rows[0].line_number > 1
    match = models.get_orbital_body_by_name(rows, "cupido")
    assert match.name == rows[0].name
    assert models.get_orbital_body_by_name(rows, "not-a-model") is None
    assert models.calc_seorbel_position(match, 2451545.0) == models.calc_orbital_position(match, 2451545.0)
    assert models.get_bundled_seorbel_path() == path


def test_all_exported_models_produce_finite_positions_and_daily_rates():
    functions = (models.calc_cupido, models.calc_hades, models.calc_zeus, models.calc_kronos,
        models.calc_apollon, models.calc_admetos, models.calc_vulkanus, models.calc_poseidon,
        models.calc_transpluto, models.calc_vulcan, models.calc_waldemath, models.calc_proserpina,
        models.calc_planet_x_pickering, models.calc_white_moon_position)
    for function in functions:
        state = function(2451545.0)
        assert isinstance(state, tuple) and len(state) == 6
        assert all(math.isfinite(value) for value in state)
        assert 0.0 <= state[0] < 360.0 and state[2] > 0.0
    assert inspect.signature(models.calc_white_moon_position).parameters["use_true_lilith"].default is False
    with pytest.raises(eph.UnknownBodyError):
        models.calc_uranian_planet(49, 2451545.0)
    with pytest.raises(ValueError):
        models.calc_hypothetical_position(-1, 2451545.0)


def test_hypothetical_models_are_deterministic_across_eight_workers():
    requests = [(body, 2451545.0 + day) for body in range(40, 48) for day in range(2)]
    expected = [models.calc_uranian_planet(*request) for request in requests]
    with ThreadPoolExecutor(max_workers=8) as pool:
        observed = list(pool.map(lambda request: models.calc_uranian_planet(*request), requests))
    assert observed == expected
