# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the private certified ellipsoid contact frame."""

from __future__ import annotations

import dataclasses
import math

import pytest

from flint import arb

from libephemeris.ellipsoid_contact import (
    _ConeSection,
    _EllipsoidContactFrame,
    _certify_regular_root_box,
    _evaluate_regular_contact,
    _regular_contact_jacobian,
)
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


def test_regular_contact_matches_exact_e1_tangency() -> None:
    """E1 satisfies ellipsoid, cone, nappe, and stationarity equations."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    evaluation = _evaluate_regular_contact(
        frame,
        cone,
        (arb(-1.6), arb(0), arb(1.2)),
        arb(1),
    )
    assert evaluation.ellipsoid.contains(0)
    assert evaluation.cone.contains(0)
    assert evaluation.nappe_radius > 0
    assert evaluation.radial_distance > 0
    assert all(component.contains(0) for component in evaluation.stationarity)


def test_regular_contact_jacobian_matches_exact_cylinder_case() -> None:
    """A dyadic cylinder gives an independently explicit 4-by-4 Jacobian."""
    frame = _frame(
        metric_km_minus_2=(
            (0.25, 0.0, 0.0),
            (0.0, 0.25, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(0.0, 0.0, 0.0),
        axis_span=(0.0, 0.0, 1.0),
    )
    jacobian = _regular_contact_jacobian(
        frame,
        _ConeSection(1.0, 1.0, 1),
        (arb(-2), arb(0), arb(0)),
        arb(1),
    )
    expected = (
        (-1.0, 0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0, -1.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 0.5, 0.0),
    )
    assert (jacobian.nrows(), jacobian.ncols()) == (4, 4)
    for row in range(4):
        for column in range(4):
            assert jacobian[row, column].contains(expected[row][column])


def test_regular_contact_jacobian_encloses_e1_finite_difference() -> None:
    """E1 analytic columns contain independent centered-difference estimates."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    point = (-1.6, 0.0, 1.2)
    jacobian = _regular_contact_jacobian(
        frame, cone, tuple(arb(value) for value in point), arb(1)
    )
    step = 1e-6
    for column in range(3):
        low = list(point)
        high = list(point)
        low[column] -= step
        high[column] += step
        lower = _evaluate_regular_contact(
            frame, cone, tuple(arb(value) for value in low), arb(1)
        )
        upper = _evaluate_regular_contact(
            frame, cone, tuple(arb(value) for value in high), arb(1)
        )
        low_values = (lower.ellipsoid,) + lower.stationarity
        high_values = (upper.ellipsoid,) + upper.stationarity
        for row, (before, after) in enumerate(zip(low_values, high_values)):
            difference = (float(after.mid()) - float(before.mid())) / (2 * step)
            assert abs(float(jacobian[row, column].mid()) - difference) < 1e-8


def test_regular_root_box_certifies_e1_kkt_root() -> None:
    """Strict Krawczyk inclusion isolates the regular E1 stationary root."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    box = (
        arb(-1.6, "0.01"),
        arb(0, "0.01"),
        arb(1.2, "0.01"),
        arb(1, "0.01"),
    )
    certificate = _certify_regular_root_box(frame, cone, box)
    assert all(
        outer.contains_interior(inner)
        for outer, inner in zip(certificate.root_box, certificate.image)
    )
    assert certificate.cone_residual.contains(0)


def test_regular_root_box_rejects_noncontracting_or_negative_lambda_box() -> None:
    """A broad image or multiplier crossing zero cannot certify a root."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    with pytest.raises(IntervalCertificationError, match="strictly inside"):
        _certify_regular_root_box(
            frame,
            cone,
            (arb(-1.6, "0.3"), arb(0, "0.3"), arb(1.2, "0.3"), arb(1.5, "0.3")),
        )
    with pytest.raises(ValueError, match="multiplier"):
        _certify_regular_root_box(
            frame,
            cone,
            (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(0, 1)),
        )


def test_regular_root_certificate_does_not_claim_cone_contact() -> None:
    """A KKT root with a changed cone radius retains a nonzero g enclosure."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(0.5, 0.8, 1)
    certificate = _certify_regular_root_box(
        frame,
        cone,
        (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(1, "0.01")),
    )
    assert not certificate.cone_residual.contains(0)


def test_regular_contact_requires_smooth_domain() -> None:
    """Axis and inactive-nappe points stay outside the regular KKT evaluator."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    with pytest.raises(IntervalCertificationError, match="radial"):
        _evaluate_regular_contact(
            frame, _ConeSection(1.0, 1.0, 1), (arb(0), arb(0), arb(0)), arb(0)
        )
    with pytest.raises(IntervalCertificationError, match="nappe"):
        _evaluate_regular_contact(
            frame,
            _ConeSection(0.0, 0.8, 1),
            (arb(1), arb(0), arb(-2)),
            arb(0),
        )


@pytest.mark.parametrize(
    "cone",
    [
        _ConeSection(-1.0, 0.8, 1),
        _ConeSection(1.0, 0.0, 1),
        _ConeSection(1.0, 1.1, 1),
        _ConeSection(1.0, 0.8, 0),
    ],
)
def test_cone_section_rejects_invalid_inputs(cone) -> None:
    """Physical cone data is never clamped or assigned a default branch."""
    with pytest.raises(ValueError):
        cone.validate()


@pytest.mark.parametrize("multiplier", [arb(-1), arb(0, 1), -1.0, None])
def test_regular_contact_rejects_uncertified_multiplier(multiplier) -> None:
    """The KKT multiplier enclosure must be wholly non-negative."""
    with pytest.raises(ValueError, match="multiplier"):
        _evaluate_regular_contact(
            _frame(),
            _ConeSection(1.0, 1.0, 1),
            (arb(-2), arb(0), arb(0)),
            multiplier,
        )


@pytest.mark.parametrize("point", [(arb(-2), 0.0, arb(0)), [arb(-2), arb(0), arb(0)]])
def test_regular_contact_rejects_malformed_point(point) -> None:
    """Wrong point containers and scalar types raise the declared ValueError."""
    with pytest.raises(ValueError, match="contact point"):
        _evaluate_regular_contact(_frame(), _ConeSection(1.0, 1.0, 1), point, arb(0))


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
