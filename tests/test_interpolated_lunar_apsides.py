"""Physical and numerical invariants for smoothed JPL lunar apsides."""

from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor

import pytest

from libephemeris.exceptions import EphemerisRangeError
from libephemeris.interpolated_lunar_apsides import (
    INTP_APOG,
    INTP_PERG,
    SUPPORTED_BODIES,
    InterpolatedLunarApsidesReader,
    _QUADRATURE_NODES,
    _QUADRATURE_WEIGHTS,
    _active_jpl_context,
    _mean_anomalistic_period_days,
    _moon_icrs_state,
    evaluate_smoothed_lunar_apsis_state,
    open_interpolated_lunar_apsides,
    state_to_polar,
)

J2000 = 2451545.0


def _unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    magnitude = math.sqrt(sum(value * value for value in vector))
    return tuple(value / magnitude for value in vector)  # type: ignore[return-value]


def test_quadrature_is_a_symmetric_unit_area_boxcar() -> None:
    assert len(_QUADRATURE_NODES) == len(_QUADRATURE_WEIGHTS) == 12
    assert sum(_QUADRATURE_WEIGHTS) == pytest.approx(2.0, abs=2e-15)
    for index in range(len(_QUADRATURE_NODES)):
        opposite = -index - 1
        assert _QUADRATURE_NODES[index] == pytest.approx(
            -_QUADRATURE_NODES[opposite], abs=2e-16
        )
        assert _QUADRATURE_WEIGHTS[index] == pytest.approx(
            _QUADRATURE_WEIGHTS[opposite], abs=2e-16
        )


@pytest.mark.parametrize("jd", [2415020.5, J2000, 2488070.5])
def test_anomalistic_period_is_derived_from_iers_mean_anomaly(jd: float) -> None:
    period = _mean_anomalistic_period_days(jd)
    assert 20.0 < period < 35.0
    assert type(period) is float


@pytest.mark.parametrize("jd", [2415020.5, J2000, 2488070.5])
def test_smoothed_apogee_and_perigee_are_antipodal(jd: float) -> None:
    with InterpolatedLunarApsidesReader(cache_size=4) as reader:
        apogee = reader.state(jd, INTP_APOG)
        perigee = reader.state(jd, INTP_PERG)

    apo_unit = _unit(apogee[:3])
    per_unit = _unit(perigee[:3])
    assert sum(a * p for a, p in zip(apo_unit, per_unit)) == pytest.approx(
        -1.0, abs=3e-15
    )
    assert math.dist(apo_unit, tuple(-value for value in per_unit)) < 3e-15
    assert math.sqrt(sum(value * value for value in apogee[:3])) > math.sqrt(
        sum(value * value for value in perigee[:3])
    )
    assert all(type(value) is float for value in apogee + perigee)


@pytest.mark.parametrize("body_id", SUPPORTED_BODIES)
def test_polar_state_is_physical_and_roundtrips_radius(body_id: int) -> None:
    with open_interpolated_lunar_apsides() as reader:
        cartesian = reader.state(J2000, body_id)
        polar = reader.polar(J2000, body_id)

    radius = math.sqrt(sum(value * value for value in cartesian[:3]))
    assert polar[2] == pytest.approx(radius, abs=1e-18)
    assert 0.0 <= polar[0] < 360.0
    assert -10.0 < polar[1] < 10.0
    assert 0.002 < polar[2] < 0.0035
    assert all(math.isfinite(value) for value in polar)


@pytest.mark.parametrize("body_id", SUPPORTED_BODIES)
def test_reported_velocity_matches_local_position_change(body_id: int) -> None:
    step = 0.01
    with InterpolatedLunarApsidesReader(cache_size=4) as reader:
        state = reader.state(J2000, body_id)
        before = reader.state(J2000 - step, body_id)
        after = reader.state(J2000 + step, body_id)

    finite_difference = tuple(
        (after[index] - before[index]) / (2.0 * step) for index in range(3)
    )
    assert state[3:] == pytest.approx(finite_difference, abs=2e-8)


def test_reader_cache_is_bounded_and_deterministic() -> None:
    reader = InterpolatedLunarApsidesReader(cache_size=2)
    first = reader.state(J2000, INTP_APOG)
    assert reader.state(J2000, INTP_APOG) is first
    reader.state(J2000 + 1.0, INTP_APOG)
    reader.state(J2000 + 2.0, INTP_APOG)
    assert len(reader._cache) == 2
    assert (J2000, INTP_APOG, False) not in reader._cache
    reader.close()


def test_reader_serializes_concurrent_access() -> None:
    reader = InterpolatedLunarApsidesReader(cache_size=4)
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(lambda _: reader.state(J2000, INTP_PERG), range(8)))
    assert all(result == results[0] for result in results)
    assert len(reader._cache) == 1
    reader.close()


def test_speed3_uses_the_same_smooth_instantaneous_state() -> None:
    with InterpolatedLunarApsidesReader() as reader:
        ordinary = reader.state(J2000, INTP_APOG)
        speed3 = reader.state(J2000, INTP_APOG, speed3=True)
    assert speed3 == ordinary


def test_injected_backend_provider_uses_the_same_definition() -> None:
    _, start_jd, end_jd, _ = _active_jpl_context()
    injected = evaluate_smoothed_lunar_apsis_state(
        J2000,
        INTP_PERG,
        _moon_icrs_state,
        jd_start=start_jd,
        jd_end=end_jd,
    )
    with InterpolatedLunarApsidesReader() as reader:
        direct = reader.polar(J2000, INTP_PERG)
    assert injected == direct


def test_reader_lifecycle_and_validation() -> None:
    reader = InterpolatedLunarApsidesReader()
    assert not reader.closed
    reader.close()
    assert reader.closed
    with pytest.raises(ValueError, match="closed"):
        reader.state(J2000, INTP_APOG)
    with pytest.raises(ValueError, match="cache_size"):
        InterpolatedLunarApsidesReader(cache_size=0)
    with pytest.raises(ValueError, match="external"):
        InterpolatedLunarApsidesReader("obsolete.bin")


def test_unsupported_body_is_rejected() -> None:
    with InterpolatedLunarApsidesReader() as reader:
        with pytest.raises(ValueError, match="Unsupported"):
            reader.state(J2000, 1)


def test_range_tracks_the_active_jpl_kernel() -> None:
    context = _active_jpl_context()
    _, start_jd, end_jd, _ = context
    with InterpolatedLunarApsidesReader() as reader:
        assert reader.jd_start == start_jd
        assert reader.jd_end == end_jd
        assert all(math.isfinite(value) for value in reader.state(start_jd, INTP_APOG))
        assert all(math.isfinite(value) for value in reader.state(end_jd, INTP_PERG))
        with pytest.raises(EphemerisRangeError):
            reader.state(math.nextafter(start_jd, -math.inf), INTP_APOG)
        with pytest.raises(EphemerisRangeError):
            reader.state(math.nextafter(end_jd, math.inf), INTP_PERG)


def test_state_to_polar_rejects_degenerate_vectors() -> None:
    with pytest.raises(Exception, match="degenerate"):
        state_to_polar((0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
