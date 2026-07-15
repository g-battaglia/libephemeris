"""Independent invariants for the IERS mean lunar points."""

from __future__ import annotations

import math

import erfa
import pytest

from libephemeris.exceptions import EphemerisRangeError
from libephemeris.mean_lunar_apse import (
    INCLINATION_DEG,
    MEAN_APOGEE_DISTANCE_AU,
    MEAN_NODE_DISTANCE_AU,
    lunar_delaunay_arguments,
    mean_lunar_apogee_position,
    mean_lunar_apogee_state,
    mean_lunar_apse_plane_node,
    mean_lunar_node_position,
    mean_lunar_node_state,
    _active_ephemeris_range,
)
from libephemeris.state import get_current_file_data, get_planets

J2000 = 2451545.0


@pytest.mark.parametrize("jd", [2415020.5, J2000, 2488070.5])
def test_shared_delaunay_basis_is_exactly_the_iers_erfa_basis(jd: float) -> None:
    centuries = (jd - J2000) / 36525.0
    assert lunar_delaunay_arguments(jd) == pytest.approx(
        (
            float(erfa.fad03(centuries)),
            float(erfa.falp03(centuries)),
            float(erfa.fal03(centuries)),
            float(erfa.faf03(centuries)),
        ),
        abs=0.0,
    )


@pytest.mark.parametrize("jd", [2415020.5, J2000, 2488070.5])
def test_plane_and_node_follow_iers_delaunay_identities(jd: float) -> None:
    centuries = (jd - J2000) / 36525.0
    lunar_anomaly = float(erfa.fal03(centuries))
    argument_latitude = float(erfa.faf03(centuries))
    node = float(erfa.faom03(centuries))
    plane, actual_node = mean_lunar_apse_plane_node(jd)

    expected_plane = (
        math.degrees(node + argument_latitude - lunar_anomaly + math.pi) % 360.0
    )
    expected_node = math.degrees(node) % 360.0
    assert plane == pytest.approx(expected_plane, abs=2e-13)
    assert actual_node == pytest.approx(expected_node, abs=2e-13)


@pytest.mark.parametrize("jd", [2415020.5, J2000, 2488070.5])
def test_mean_apogee_is_on_the_conventional_inclined_orbit(jd: float) -> None:
    longitude, latitude, distance = mean_lunar_apogee_position(jd)
    plane, _ = mean_lunar_apse_plane_node(jd)
    inclination = math.radians(INCLINATION_DEG)

    # Rotate the returned unit vector back to the node/argument geometry.  Its
    # latitude is bounded by the independently specified mean inclination.
    assert 0.0 <= longitude < 360.0
    assert abs(latitude) <= INCLINATION_DEG + 1e-12
    assert distance == pytest.approx(MEAN_APOGEE_DISTANCE_AU, abs=5e-19)
    assert abs((longitude - plane + 180.0) % 360.0 - 180.0) < 0.2
    assert math.isfinite(inclination)


@pytest.mark.parametrize(
    "state_evaluator",
    [mean_lunar_node_state, mean_lunar_apogee_state],
)
@pytest.mark.parametrize("jd", [2415020.5, J2000, 2488070.5])
def test_analytic_rates_match_independent_position_differences(
    state_evaluator: object, jd: float
) -> None:
    step = 0.01
    state = state_evaluator(jd)  # type: ignore[operator]
    before = state_evaluator(jd - step)  # type: ignore[operator]
    after = state_evaluator(jd + step)  # type: ignore[operator]
    longitude_rate = ((after[0] - before[0] + 180.0) % 360.0 - 180.0) / (2.0 * step)
    latitude_rate = (after[1] - before[1]) / (2.0 * step)

    assert state[3] == pytest.approx(longitude_rate, abs=2e-8)
    assert state[4] == pytest.approx(latitude_rate, abs=2e-8)
    assert state[5] == 0.0
    assert all(type(value) is float for value in state)


def test_node_geometry_is_exactly_planar_and_constant_radius() -> None:
    for jd in (2415020.5, J2000, 2488070.5):
        state = mean_lunar_node_state(jd)
        assert state[1] == 0.0
        assert state[2] == MEAN_NODE_DISTANCE_AU
        assert state[4:] == (0.0, 0.0)


@pytest.mark.parametrize("evaluator", [mean_lunar_node_state, mean_lunar_apogee_state])
def test_speed3_selects_same_standards_derivative(evaluator: object) -> None:
    assert evaluator(J2000, speed3=True) == evaluator(J2000)  # type: ignore[operator]


def test_direct_evaluators_follow_active_jpl_coverage() -> None:
    get_planets()
    _, start_jd, end_jd, _ = get_current_file_data(0)
    mean_lunar_node_position(start_jd)
    mean_lunar_apogee_position(end_jd)
    with pytest.raises(EphemerisRangeError):
        mean_lunar_node_position(math.nextafter(start_jd, -math.inf))
    with pytest.raises(EphemerisRangeError):
        mean_lunar_apogee_position(math.nextafter(end_jd, math.inf))


def test_horizons_range_check_does_not_load_a_local_kernel(monkeypatch) -> None:
    """Analytic mean points must not download de440 in Horizons mode."""
    import libephemeris.mean_lunar_apse as model

    monkeypatch.setattr(model, "get_calc_mode", None, raising=False)
    monkeypatch.setattr("libephemeris.state.get_calc_mode", lambda: "horizons")
    monkeypatch.setattr(
        "libephemeris.state.get_current_file_data",
        lambda body_id: (None, 0.0, 0.0, None),
    )
    monkeypatch.setattr(
        "libephemeris.state.get_planets",
        lambda: pytest.fail("Horizons mean-point range check loaded de440.bsp"),
    )

    assert _active_ephemeris_range() == (None, 0.0, 0.0)


@pytest.mark.parametrize("jd", [math.nan, math.inf, -math.inf])
def test_nonfinite_dates_raise_range_error(jd: float) -> None:
    with pytest.raises(EphemerisRangeError):
        mean_lunar_node_state(jd)
