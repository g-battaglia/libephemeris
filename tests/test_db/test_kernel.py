# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Direct tests of the portable kernel: records in, primitive results out."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import libephemeris as ephemeris
from libephemeris.db.contract import Series, decode_segment
from libephemeris.db.kernel import (
    evaluate_body,
    evaluate_nutation,
    interpolate_delta_t,
    required_segment_keys,
    segment_coordinates,
)
from libephemeris.exceptions import DBDataError


@pytest.mark.parametrize(
    "jd, position, velocity", [(10.0, 2.0, -10.0), (11.0, -2.0, 2.0), (12.0, 6.0, 14.0)]
)
def test_portable_fixture_evaluates_without_an_operation(jd, position, velocity):
    """Exercise explicit functions without a store, owner or current calc mode.

    Args:
        jd: Synthetic epoch, including both segment endpoints.
        position: Exact first polynomial value from 6*tau**2 + 2*tau - 2.
        velocity: Exact first derivative with the fixture's interval scale.
    """
    fixture_path = Path(__file__).parent / "fixtures/contract_v1.json"
    fixture = json.loads(fixture_path.read_text())
    series = Series(**fixture["series"])
    coefficients = decode_segment(series, bytes.fromhex(fixture["payload_hex"]))
    index, tau = segment_coordinates(series, jd)
    assert index == 0
    result = evaluate_body(series, coefficients, tau)
    assert result[0][0] == position
    assert result[1][0] == velocity
    assert result[0][2] == 7.0
    assert result[1][2] == 0.0
    assert all(type(value) is float for vector in result for value in vector)


def test_spherical_longitude_uses_euclidean_modulo():
    """Document the difference between signed remainder and Euclidean modulo."""
    series = Series(0, 1, 1, 10.0, 12.0, 2.0, 0, 3)
    position, velocity = evaluate_body(series, (-2.0, 0.0, 1.0), 0.0)
    assert position == (358.0, 0.0, 1.0)
    assert velocity == (0.0, 0.0, 0.0)


def test_nutation_is_an_explicit_two_component_evaluation():
    """Reuse the native recurrence with synthetic component-major inputs."""
    series = Series(-1, 0, 1, 10.0, 14.0, 4.0, 2, 2)
    assert evaluate_nutation(series, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0), 0.0) == (
        -2.0,
        -2.0,
    )


@pytest.mark.parametrize(
    "flags, dependencies",
    [
        (
            0,
            [
                -1,
                ephemeris.SUN,
                ephemeris.MARS,
                ephemeris.JUPITER,
                ephemeris.SATURN,
                ephemeris.EARTH,
            ],
        ),
        (ephemeris.FLG_BARYCTR | ephemeris.FLG_NONUT, [ephemeris.MARS]),
        (ephemeris.FLG_HELCTR | ephemeris.FLG_J2000, [ephemeris.SUN, ephemeris.MARS]),
    ],
)
def test_prefetch_plan_is_explicit_and_deterministic(flags, dependencies):
    """Compare key planning independently of PostgreSQL and mutable buffers.

    Args:
        flags: Normalized calculation flags.
        dependencies: Sorted body identifiers expected in the plan.
    """
    inventory = {}
    for body in [
        -1,
        ephemeris.SUN,
        ephemeris.MARS,
        ephemeris.JUPITER,
        ephemeris.SATURN,
        ephemeris.EARTH,
    ]:
        inventory[body] = Series(body, 0, 2, 10.0, 14.0, 2.0, 2, 2 if body == -1 else 3)
    expected = [(body, index) for body in dependencies for index in (0, 1)]
    assert required_segment_keys(inventory, 12.0, ephemeris.MARS, flags) == expected
    assert required_segment_keys(inventory, 20.0, ephemeris.MARS, flags) == []


def test_delta_t_interpolation_and_missing_optional_samples():
    """Retain interpolation/clamping while rejecting malformed persisted data."""
    assert interpolate_delta_t([(10.0, 0.1), (14.0, 0.2)], 12.0) == pytest.approx(0.15)
    assert interpolate_delta_t([(10.0, 0.1)], 9.0) == 0.1
    with pytest.raises(ValueError, match="No Delta-T"):
        interpolate_delta_t([], 12.0)
    with pytest.raises(DBDataError, match="ordered"):
        interpolate_delta_t([(10.0, 0.1), (10.0, 0.2)], 12.0)
    with pytest.raises(DBDataError, match="samples"):
        interpolate_delta_t([(10.0, math.nan)], 12.0)
