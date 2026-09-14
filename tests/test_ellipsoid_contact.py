# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the private certified ellipsoid contact frame."""

from __future__ import annotations

import dataclasses
import math

import pytest

from libephemeris.ellipsoid_contact import _EllipsoidContactFrame
from libephemeris.intervals import IntervalCertificationError


def _frame(**overrides) -> _EllipsoidContactFrame:
    values = {
        "metric_km_minus_2": (
            (0.25, 0.0, 0.0),
            (0.0, 0.25, 0.0),
            (0.0, 0.0, 0.25),
        ),
        "axis_point_km": (-3.5, 2.0, 6.0),
        "axis_span": (0.0, 3.0, 4.0),
        "penumbra_branch_sign": 1,
        "core_branch_sign": -1,
    }
    values.update(overrides)
    return _EllipsoidContactFrame(**values)


def test_frame_is_private_frozen_and_equal_by_value() -> None:
    """The internal frame is owned immutable data, not ambient metadata."""
    frame = _frame()
    assert frame == _frame()
    with pytest.raises(dataclasses.FrozenInstanceError):
        frame.axis_span = (1.0, 0.0, 0.0)  # type: ignore[misc]


def test_axis_line_derives_unit_direction_and_closest_anchor() -> None:
    """Normalization and orthogonality are construction identities in Arb."""
    line = _frame().certified_axis_line()
    norm_squared = sum((value * value for value in line.direction), 0)
    dot = sum(
        (line.direction[index] * line.anchor_km[index] for index in range(3)),
        0,
    )
    assert norm_squared.contains(1)
    assert dot.contains(0)
    assert line.direction[0].is_zero()
    assert line.direction[1].contains(0.6)
    assert line.direction[2].contains(0.8)
    assert all(value.is_finite() for value in line.anchor_km)


def test_axis_line_is_invariant_under_another_point_on_the_line() -> None:
    """Changing the supplied point along the same line preserves its anchor."""
    first = _frame(axis_point_km=(-3.5, 2.0, 6.0)).certified_axis_line()
    second = _frame(axis_point_km=(-3.5, 5.0, 10.0)).certified_axis_line()
    for left, right in zip(first.anchor_km, second.anchor_km):
        assert left.overlaps(right)


def test_metric_matrix_preserves_exact_binary64_inputs() -> None:
    """The metric enters Arb without decimal re-rounding."""
    metric = _frame().metric_ball_matrix()
    assert metric[0, 0].is_exact()
    assert metric[0, 0] == 0.25
    assert metric[2, 2] == 0.25


@pytest.mark.parametrize("sign", [0, 2, -2, True, 1.0])
def test_branch_signs_are_exactly_plus_or_minus_one(sign) -> None:
    """Nappe selection admits no Boolean, float, zero, or wider integer."""
    with pytest.raises(ValueError, match="sign"):
        _frame(core_branch_sign=sign).validate()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("axis_span", (0.0, 0.0, 0.0), "squared norm"),
        ("axis_span", (math.inf, 0.0, 1.0), "finite"),
        ("axis_point_km", (0.0, math.nan, 0.0), "finite"),
        ("axis_point_km", [0.0, 0.0, 0.0], "tuple"),
    ],
)
def test_axis_vectors_fail_closed(field, value, message) -> None:
    """Malformed or non-finite source vectors never enter Arb geometry."""
    with pytest.raises((ValueError, IntervalCertificationError), match=message):
        _frame(**{field: value}).validate()


def test_metric_rejects_asymmetry_and_nonpositive_definiteness() -> None:
    """The frame does not symmetrize or regularize an invalid ellipsoid."""
    with pytest.raises(ValueError, match="symmetric"):
        _frame(
            metric_km_minus_2=(
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        ).validate()
    with pytest.raises(IntervalCertificationError, match="pivot"):
        _frame(
            metric_km_minus_2=(
                (1.0, 0.0, 0.0),
                (0.0, -1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        ).validate()
