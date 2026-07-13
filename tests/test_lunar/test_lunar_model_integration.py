"""Runtime integration tests for standards/JPL lunar points."""

from __future__ import annotations

from collections.abc import Callable, Iterator

import pytest

import libephemeris as swe
from libephemeris import lunar
from libephemeris.interpolated_lunar_apsides import INTP_APOG, INTP_PERG

J2000 = 2451545.0


@pytest.fixture
def skyfield_mode() -> Iterator[None]:
    previous = swe.get_calc_mode()
    swe.set_calc_mode("skyfield")
    try:
        yield
    finally:
        swe.close()
        swe.set_calc_mode(previous)


def test_mean_public_helpers_are_thin_standards_wrappers() -> None:
    node_state = lunar.calc_mean_lunar_node_state(J2000)
    apogee_state = lunar.calc_mean_lilith_state(J2000)

    assert lunar.calc_mean_lunar_node(J2000) == node_state[0]
    assert lunar.calc_mean_lilith(J2000) == apogee_state[0]
    assert lunar.calc_mean_lilith_with_latitude(J2000) == apogee_state[:2]
    assert node_state[1] == 0.0
    assert node_state[4:] == (0.0, 0.0)
    assert apogee_state[5] == 0.0
    assert all(type(value) is float for value in node_state + apogee_state)


@pytest.mark.parametrize(
    ("body_id", "position_helper"),
    [
        (INTP_APOG, lunar.calc_interpolated_apogee),
        (INTP_PERG, lunar.calc_interpolated_perigee),
    ],
)
def test_interpolated_helpers_share_the_jpl_smoothed_state(
    body_id: int,
    position_helper: Callable[[float], tuple[float, float, float]],
) -> None:
    state = lunar.calc_interpolated_apse_state(J2000, body_id)
    assert position_helper(J2000) == state[:3]
    assert 0.0 <= state[0] < 360.0
    assert -10.0 < state[1] < 10.0
    assert 0.002 < state[2] < 0.0035
    assert all(type(value) is float for value in state)


def test_interpolated_state_rejects_non_apse_body() -> None:
    with pytest.raises(ValueError, match="Unsupported"):
        lunar.calc_interpolated_apse_state(J2000, swe.MOON)


def test_evaluator_lifecycle_distinguishes_reset_from_close() -> None:
    lunar.release_interpolated_apsides_reader()
    evaluator = lunar._get_interpolated_apsides_reader()
    evaluator.polar(J2000, INTP_APOG)

    swe.reset_session()
    assert lunar._get_interpolated_apsides_reader() is evaluator
    assert not evaluator.closed

    swe.close()
    assert evaluator.closed
    reopened = lunar._get_interpolated_apsides_reader()
    assert reopened is not evaluator
    assert not reopened.closed
    lunar.release_interpolated_apsides_reader()


@pytest.mark.parametrize(
    ("body_id", "state_helper"),
    [
        (swe.MEAN_NODE, lunar.calc_mean_lunar_node_state),
        (swe.MEAN_APOG, lunar.calc_mean_lilith_state),
    ],
)
def test_calc_nonut_speed_uses_the_mean_standards_state(
    body_id: int,
    state_helper: Callable[[float], tuple[float, float, float, float, float, float]],
    skyfield_mode: None,
) -> None:
    result, _ = swe.calc(J2000, body_id, swe.FLG_NONUT | swe.FLG_SPEED)
    assert result == pytest.approx(state_helper(J2000), abs=2e-12)
    assert all(type(value) is float for value in result)


@pytest.mark.parametrize("body_id", [swe.INTP_APOG, swe.INTP_PERG])
def test_calc_speed_uses_the_smoothed_jpl_state(
    body_id: int, skyfield_mode: None
) -> None:
    result, _ = swe.calc(J2000, body_id, swe.FLG_SPEED)
    expected = lunar.calc_interpolated_apse_state(J2000, body_id)
    assert result == pytest.approx(expected, abs=2e-12)
    assert all(type(value) is float for value in result)


@pytest.mark.parametrize(
    "body_id", [swe.MEAN_NODE, swe.MEAN_APOG, swe.INTP_APOG, swe.INTP_PERG]
)
def test_calc_without_speed_zeroes_velocity_channels(
    body_id: int, skyfield_mode: None
) -> None:
    result, _ = swe.calc(J2000, body_id, swe.FLG_NONUT)
    assert result[3:] == (0.0, 0.0, 0.0)
