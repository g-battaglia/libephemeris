# SPDX-License-Identifier: AGPL-3.0-only
"""Synthetic tests for lunar-eclipse maximum refinement."""

from __future__ import annotations

import math
from collections.abc import Callable

import pytest

from libephemeris import eclipse
from libephemeris.constants import FLG_EQUATORIAL, FLG_MOSEPH, FLG_SPEED, FLG_XYZ
from libephemeris.constants import MOON, SUN
from libephemeris.exceptions import ConvergenceError, EphemerisRangeError


def _install_synthetic_states(
    monkeypatch: pytest.MonkeyPatch,
    seed: float,
    angle: Callable[[float], float],
    sun_distance: Callable[[float], float] = lambda _t: 1.0,
) -> list[tuple[float, int, int]]:
    """Installa stati apparenti coerenti e restituisce il registro calls."""
    calls: list[tuple[float, int, int]] = []
    moon_distance = 0.00257

    def fake_calc_ut(t: float, body: int, flags: int):
        calls.append((t, body, flags))
        if body == MOON:
            position = (moon_distance, 0.0, 0.0)
        elif body == SUN:
            alpha = angle(t - seed)
            radius = sun_distance(t - seed)
            position = (
                moon_distance - radius * math.cos(alpha),
                radius * math.sin(alpha),
                0.0,
            )
        else:  # pragma: no cover - make the two-state contract explicit
            raise AssertionError(f"unexpected body: {body}")
        return ((*position, 0.0, 0.0, 0.0), flags)

    monkeypatch.setattr(eclipse, "calc_ut", fake_calc_ut)
    return calls


def test_maximizes_apparent_overlap_within_seed_window(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A varying solar radius moves the maximum away from minimum separation."""
    seed = 2460000.0
    curvature = 0.1
    distance_slope = -0.5
    calls = _install_synthetic_states(
        monkeypatch,
        seed,
        lambda dt: curvature * dt * dt,
        lambda dt: 1.0 + distance_slope * dt,
    )

    actual = eclipse._lun_eclipse_max_time(seed, FLG_MOSEPH | FLG_SPEED)

    # The published objective must rise toward the result and fall after it;
    # this also demonstrates displacement from minimum separation at the seed.
    assert isinstance(actual, float)
    assert actual > seed + 0.01
    values = []
    for offset in (-1.0e-5, 0.0, 1.0e-5):
        dt = actual + offset - seed
        distance = 1.0 + distance_slope * dt
        values.append(math.asin(eclipse._ECL_RSUN_AU / distance) - curvature * dt * dt)
    assert values[1] > values[0]
    assert values[1] > values[2]
    assert calls
    assert all(seed - 0.3 <= t <= seed + 0.3 for t, _corpo, _flags in calls)
    assert {body for _t, body, _flags in calls} == {MOON, SUN}
    state_flags = FLG_MOSEPH | FLG_EQUATORIAL | FLG_XYZ
    assert all(flags == state_flags for _t, _corpo, flags in calls)


@pytest.mark.parametrize("offset", [-0.3, 0.3])
def test_includes_both_closed_window_endpoints(
    monkeypatch: pytest.MonkeyPatch, offset: float
) -> None:
    """A unique maximum at either endpoint is returned exactly."""
    seed = 2460100.0
    endpoint = seed + offset
    _install_synthetic_states(
        monkeypatch,
        seed,
        lambda dt: 0.4 * abs((seed + dt) - endpoint),
    )

    actual = eclipse._lun_eclipse_max_time(seed)

    assert actual == endpoint


@pytest.mark.parametrize("target", [-0.295, 0.295])
@pytest.mark.parametrize("exponent", [1, 2])
def test_finds_maximum_in_segment_adjacent_to_endpoint(
    monkeypatch: pytest.MonkeyPatch, target: float, exponent: int
) -> None:
    """The first and last segments are not mistaken for the endpoints."""
    seed = 2460000.0
    _install_synthetic_states(
        monkeypatch,
        seed,
        lambda dt: 0.4 * abs(dt - target) ** exponent,
    )

    actual = eclipse._lun_eclipse_max_time(seed)

    assert actual == pytest.approx(seed + target, abs=2.0e-8)


@pytest.mark.parametrize(
    "target",
    [
        sign * (0.275 + index * 1.0e-9)
        for sign in (-1.0, 1.0)
        for index in range(-20, 21)
    ],
)
def test_does_not_duplicate_unique_peak_at_boundary_transition(
    monkeypatch: pytest.MonkeyPatch, target: float
) -> None:
    """Interior and boundary brackets describe one candidate region."""
    seed = 2460000.0
    _install_synthetic_states(
        monkeypatch,
        seed,
        lambda dt: math.hypot(1.0e-5, 0.01 * (dt - target)),
    )

    actual = eclipse._lun_eclipse_max_time(seed)

    assert actual == pytest.approx(seed + target, abs=2.0e-9)


def test_rejects_two_separate_maxima(monkeypatch: pytest.MonkeyPatch) -> None:
    """A sampled valley preserves two independent maximum regions."""
    seed = 2460000.0
    _install_synthetic_states(
        monkeypatch,
        seed,
        lambda dt: 0.2 * min(abs(dt + 0.1), abs(dt - 0.1)),
    )

    with pytest.raises(ConvergenceError, match="not unique"):
        eclipse._lun_eclipse_max_time(seed)


def test_propagates_missing_coverage_at_endpoint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Required coverage spans the complete window, endpoints included."""
    seed = 2460200.0
    limit = seed - 0.3
    error = EphemerisRangeError("synthetic coverage unavailable", requested_jd=limit)

    def fake_calc_ut(t: float, body: int, flags: int):
        if t == limit and body == SUN:
            raise error
        if body == MOON:
            position = (0.00257, 0.0, 0.0)
        else:
            position = (-0.99743, t - seed, 0.0)
        return ((*position, 0.0, 0.0, 0.0), flags)

    monkeypatch.setattr(eclipse, "calc_ut", fake_calc_ut)

    with pytest.raises(EphemerisRangeError) as caught:
        eclipse._lun_eclipse_max_time(seed)

    assert caught.value is error


def test_rejects_nonunique_maximum_and_nonfinite_seeds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A constant objective does not permit choosing an arbitrary instant."""
    seed = 2460300.0
    _install_synthetic_states(monkeypatch, seed, lambda _dt: 0.0)

    with pytest.raises(ConvergenceError, match="not unique"):
        eclipse._lun_eclipse_max_time(seed)

    for nonfinite_seed in (math.nan, math.inf, -math.inf):
        with pytest.raises(ValueError, match="finite"):
            eclipse._lun_eclipse_max_time(nonfinite_seed)
