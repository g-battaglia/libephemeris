# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the private certified ellipsoid contact frame."""

from __future__ import annotations

import dataclasses
import math
from fractions import Fraction

import pytest

from flint import arb, ctx

from libephemeris.ellipsoid_contact import (
    _ConeSection,
    _CrossDualWitness,
    _EllipsoidContactFrame,
    _GlobalReachStatus,
    _ProjectedDualWitness,
    _RegularRootCertificate,
    _certify_regular_root_box,
    _evaluate_cone_apex,
    _evaluate_nappe_boundary,
    _evaluate_radial_subgradient,
    _evaluate_global_reach,
    _evaluate_regular_contact,
    _nappe_boundary_jacobian,
    _reduce_regular_residuals,
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


def _dual(
    span=(Fraction(0), Fraction(0), Fraction(1)),
    seed=(Fraction(1), Fraction(0), Fraction(0)),
    scale=Fraction(1),
    multiplier=Fraction(0),
) -> _ProjectedDualWitness:
    return _ProjectedDualWitness(span, seed, scale, multiplier)


def test_global_reach_keeps_unproved_e1_equality_unresolved() -> None:
    """Matching zero-containing bounds are not an exact contact proof."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    proof = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        (arb(-1.6), arb(0), arb(1.2)),
        _dual(scale=Fraction(4, 5)),
    )
    assert proof.status is _GlobalReachStatus.UNRESOLVED


def test_global_reach_proves_e4_miss() -> None:
    """The exact E4 support witness gives the global lower bound +1."""
    root_two = math.sqrt(2.0)
    frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(-5.0, 0.0, 0.0),
        axis_span=(0.0, root_two, root_two),
    )
    span = tuple(Fraction.from_float(value) for value in frame.axis_span)
    proof = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 1.0, 1),
        None,
        _ProjectedDualWitness(
            span,
            (Fraction(1), Fraction(0), Fraction(0)),
            Fraction(1) / _dot_fraction(span, span),
            Fraction(0),
        ),
    )
    assert proof.status is _GlobalReachStatus.MISS
    assert proof.lower_bound > 0  # type: ignore[operator]


def _dot_fraction(left, right) -> Fraction:
    return sum((a * b for a, b in zip(left, right)), Fraction(0))


def test_global_reach_uses_strict_primal_upper_bound() -> None:
    """A feasible point strictly inside a zero-radius cone proves reach."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    proof = _evaluate_global_reach(
        frame,
        _ConeSection(0.0, 0.8, 1),
        (arb(0), arb(0), arb(1)),
        None,
    )
    assert proof.status is _GlobalReachStatus.REACH
    assert proof.upper_bound < 0  # type: ignore[operator]


def test_global_reach_classifies_internal_certification_failure_unresolved(
    monkeypatch,
) -> None:
    """A proof-engine failure is not mislabeled as invalid source data."""
    monkeypatch.setattr(
        "libephemeris.ellipsoid_contact._dual_vector",
        lambda *_: (_ for _ in ()).throw(IntervalCertificationError("synthetic")),
    )
    proof = _evaluate_global_reach(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 0.8, 1),
        None,
        _dual(),
    )
    assert proof.status is _GlobalReachStatus.UNRESOLVED
    assert proof.reason == "synthetic"


def test_cross_dual_payload_proves_same_e4_miss() -> None:
    """The exact cross-product witness implements the second approved variant."""
    frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=(-5.0, 0.0, 0.0),
        axis_span=(0.0, 1.0, 1.0),
    )
    span = tuple(Fraction.from_float(value) for value in frame.axis_span)
    witness = _CrossDualWitness(
        span,
        (Fraction(0), Fraction(0), Fraction(-1)),
        Fraction(-1),
        Fraction(0),
    )
    proof = _evaluate_global_reach(frame, _ConeSection(1.0, 1.0, 1), None, witness)
    assert proof.status is _GlobalReachStatus.MISS
    assert proof.lower_bound > 0  # type: ignore[operator]


def test_global_reach_rejects_nonpoint_primal_and_missing_witnesses() -> None:
    """A primal witness is an exact point and at least one witness is required."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    nonpoint = _evaluate_global_reach(
        frame,
        _ConeSection(0.0, 0.8, 1),
        (arb(0, 1), arb(0), arb(1)),
        None,
    )
    assert nonpoint.status is _GlobalReachStatus.INVALID
    missing = _evaluate_global_reach(frame, _ConeSection(0.0, 0.8, 1), None, None)
    assert missing.status is _GlobalReachStatus.INVALID


def test_global_reach_rejects_wrong_dual_span_or_norm() -> None:
    """The exact dual payload stays bound to the frame and norm radius."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    wrong_span = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        (arb(0), arb(0), arb(0)),
        _dual(span=(Fraction(0), Fraction(1), Fraction(0))),
    )
    assert wrong_span.status is _GlobalReachStatus.INVALID
    too_large = _evaluate_global_reach(
        frame,
        _ConeSection(1.0, 0.8, 1),
        None,
        _dual(scale=Fraction(2)),
    )
    assert too_large.status is _GlobalReachStatus.UNRESOLVED


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


@pytest.mark.parametrize(
    ("axis_point", "axis_span", "expected_contact"),
    [
        ((-4.0, 0.0, 0.0), (0.0, 0.0, 1.0), True),
        ((-4.0, 0.0, 0.0), (0.0, 1.0, 1.0), True),
        ((-5.0, 0.0, 0.0), (0.0, 1.0, 1.0), False),
    ],
)
def test_regular_root_box_covers_e2_e3_and_e4(
    axis_point, axis_span, expected_contact
) -> None:
    """Principal/inclined tangencies and the inclined miss stay distinct."""
    frame = _frame(
        metric_km_minus_2=(
            (1.0 / 9.0, 0.0, 0.0),
            (0.0, 1.0 / 9.0, 0.0),
            (0.0, 0.0, 0.25),
        ),
        axis_point_km=axis_point,
        axis_span=axis_span,
    )
    certificate = _certify_regular_root_box(
        frame,
        _ConeSection(1.0, 1.0, 1),
        (arb(-3, "0.001"), arb(0, "0.001"), arb(0, "0.001"), arb(1.5, "0.001")),
    )
    if expected_contact:
        assert certificate.cone_residual.contains(0)
    else:
        assert certificate.cone_residual > 0


def test_regular_residual_reduction_orders_supplied_candidates_only() -> None:
    """The reducer selects a separated minimum without claiming completeness."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    first = _certify_regular_root_box(
        frame,
        cone,
        (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(1, "0.01")),
    )
    second = _RegularRootCertificate(
        frame,
        cone,
        first.root_box,
        first.image,
        arb(1.0, "0.01"),
    )
    with pytest.raises(IntervalCertificationError, match="sign"):
        _reduce_regular_residuals(cone, (first, second))
    separated = _RegularRootCertificate(
        frame,
        cone,
        first.root_box,
        first.image,
        arb(-0.5, "0.01"),
    )
    reduction = _reduce_regular_residuals(cone, (separated, second))
    assert reduction.minimum is separated.cone_residual
    assert reduction.candidate_count == 2
    assert reduction.sign == 1


def test_regular_residual_reduction_rejects_overlapping_candidates() -> None:
    """Overlapping residual ranges cannot establish one ordered candidate."""
    frame = _frame(axis_point_km=(-3.5, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    cone = _ConeSection(1.0, 0.8, 1)
    first = _certify_regular_root_box(
        frame,
        cone,
        (arb(-1.6, "0.01"), arb(0, "0.01"), arb(1.2, "0.01"), arb(1, "0.01")),
    )
    with pytest.raises(IntervalCertificationError, match="minimum"):
        _reduce_regular_residuals(cone, (first, first))
    with pytest.raises(ValueError, match="non-empty"):
        _reduce_regular_residuals(cone, ())
    other_cone = _ConeSection(0.5, 0.8, 1)
    with pytest.raises(ValueError, match="one frame and cone"):
        _reduce_regular_residuals(other_cone, (first,))


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


def test_nappe_boundary_evaluator_accepts_zero_containing_radius() -> None:
    """The active r=0 equation is returned instead of forced into r>0."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    evaluation = _evaluate_nappe_boundary(
        frame,
        _ConeSection(1.0, 0.8, 1),
        (arb(1), arb(0), arb(-4) / 3),
        arb(0),
        arb(0),
    )
    assert evaluation.nappe_radius.contains(0)
    assert evaluation.radial_distance > 0
    assert all(component.is_finite() for component in evaluation.stationarity)


def test_nappe_boundary_jacobian_matches_centered_cylinder() -> None:
    """A zero-angle active plane gives an independently explicit Jacobian."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    jacobian = _nappe_boundary_jacobian(
        frame,
        _ConeSection(0.0, 1.0, 1),
        (arb(-2), arb(0), arb(0)),
        arb(1),
        arb(0),
    )
    expected = (
        (-1.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0),
        (0.5, 0.0, 0.0, -1.0, 0.0),
        (0.0, 1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.5, 0.0, 0.0),
    )
    assert (jacobian.nrows(), jacobian.ncols()) == (5, 5)
    for row in range(5):
        for column in range(5):
            assert jacobian[row, column].contains(expected[row][column])


def test_nappe_boundary_rejects_uncertified_multiplier_or_radial_axis() -> None:
    """The active set keeps non-negative mu and smooth rho contracts."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    with pytest.raises(ValueError, match="nappe multiplier"):
        _evaluate_nappe_boundary(
            frame,
            _ConeSection(1.0, 0.8, 1),
            (arb(1), arb(0), arb(-4) / 3),
            arb(0),
            arb(0, 1),
        )
    with pytest.raises(IntervalCertificationError, match="radial"):
        _evaluate_nappe_boundary(
            frame,
            _ConeSection(0.0, 0.8, 1),
            (arb(0), arb(0), arb(0)),
            arb(0),
            arb(0),
        )


def test_cone_apex_solves_joint_radial_and_nappe_equations() -> None:
    """The nonzero-angle apex is derived directly without a radial offset."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    apex = _evaluate_cone_apex(frame, _ConeSection(1.0, 0.8, 1))
    assert apex.axial_parameter.contains(-4.0 / 3.0)
    assert apex.point_km[0].is_zero()
    assert apex.point_km[1].is_zero()
    assert apex.point_km[2].contains(-4.0 / 3.0)
    assert apex.cone.contains(0)
    assert apex.ellipsoid < 0


def test_cone_apex_is_finite_in_offset_tilted_frame() -> None:
    """Construction identities avoid dependent interval subtraction at the axis."""
    apex = _evaluate_cone_apex(_frame(), _ConeSection(1.0, 0.8, -1))
    assert all(value.is_finite() for value in apex.point_km)
    assert apex.ellipsoid.is_finite()
    assert apex.cone.is_exact()
    assert apex.cone.is_zero()


def test_cone_apex_rejects_zero_or_unresolved_angle() -> None:
    """A cylinder has no finite apex and unresolved sine fails closed."""
    with pytest.raises(ValueError, match="no finite apex"):
        _evaluate_cone_apex(_frame(), _ConeSection(1.0, 1.0, 1))
    previous = ctx.prec
    try:
        ctx.prec = 32
        with pytest.raises(IntervalCertificationError, match="sine"):
            _evaluate_cone_apex(
                _frame(), _ConeSection(1.0, math.nextafter(1.0, 0.0), 1)
            )
    finally:
        ctx.prec = previous


def test_radial_subgradient_evaluates_axis_witness_constraints() -> None:
    """A centered point and zero witness expose exact nonsmooth residuals."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    evaluation = _evaluate_radial_subgradient(
        frame,
        _ConeSection(1.0, 1.0, 1),
        (arb(0), arb(0), arb(0)),
        (arb(0), arb(0), arb(0)),
        arb(0),
        arb(0),
    )
    assert evaluation.radial_squared.is_zero()
    assert evaluation.witness_axis_dot.is_zero()
    assert evaluation.witness_norm_squared.is_zero()
    assert evaluation.nappe_radius == 1
    assert evaluation.ellipsoid == -1
    assert all(component.is_zero() for component in evaluation.stationarity)


def test_radial_subgradient_retains_witness_inequality_residual() -> None:
    """An inadmissible witness is reported, not silently normalized or clamped."""
    evaluation = _evaluate_radial_subgradient(
        _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0)),
        _ConeSection(1.0, 1.0, 1),
        (arb(0), arb(0), arb(0)),
        (arb(2), arb(0), arb(0)),
        arb(0),
        arb(0),
    )
    assert evaluation.witness_norm_squared == 4
    assert not evaluation.witness_norm_squared <= 1


def test_radial_subgradient_rejects_bad_witness_or_multiplier() -> None:
    """Nonsmooth KKT inputs retain exact Arb and non-negative contracts."""
    frame = _frame(axis_point_km=(0.0, 0.0, 0.0), axis_span=(0.0, 0.0, 1.0))
    with pytest.raises(ValueError, match="witness"):
        _evaluate_radial_subgradient(
            frame,
            _ConeSection(1.0, 1.0, 1),
            (arb(0), arb(0), arb(0)),
            (arb(0), 0.0, arb(0)),
            arb(0),
            arb(0),
        )
    with pytest.raises(ValueError, match="nappe multiplier"):
        _evaluate_radial_subgradient(
            frame,
            _ConeSection(1.0, 1.0, 1),
            (arb(0), arb(0), arb(0)),
            (arb(0), arb(0), arb(0)),
            arb(0),
            arb(0, 1),
        )


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
