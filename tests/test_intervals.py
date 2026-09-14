# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the certified Arb interval primitives."""

from __future__ import annotations

import math

import pytest
from flint import arb, arb_mat, ctx

from libephemeris.intervals import (
    IntervalCertificationError,
    ball_from_bounds,
    ball_from_float,
    best_root_float,
    certified_float,
    certified_sign,
    interval_cholesky,
    interval_precision,
    isolate_unique_root,
    krawczyk_image,
    strictly_contains,
    strictly_contains_vector,
)


@pytest.mark.parametrize(
    "value",
    [0.0, -0.0, 0.1, -1.25, math.pi, math.nextafter(1.0, math.inf)],
)
def test_ball_from_float_contains_exact_binary64(value: float) -> None:
    """Binary64 inputs enter Arb without a decimal round trip."""
    ball = ball_from_float(value)
    assert ball.is_exact()
    assert float(ball) == value


def test_ball_from_bounds_encloses_both_endpoints() -> None:
    """A bound constructor includes both exact endpoint values."""
    interval = ball_from_bounds(0.1, 0.2)
    assert interval.contains(ball_from_float(0.1))
    assert interval.contains(ball_from_float(0.2))


def test_precision_context_restores_the_previous_setting() -> None:
    """Temporary precision changes do not leak into later calculations."""
    previous = ctx.prec
    with interval_precision(192):
        assert ctx.prec == 192
    assert ctx.prec == previous


@pytest.mark.parametrize(
    ("value", "expected"), [(arb(2), 1), (arb(-2), -1), (arb(0), 0)]
)
def test_certified_sign_accepts_separated_or_exact_zero(
    value: arb, expected: int
) -> None:
    """Separated balls and exact zero have deterministic signs."""
    assert certified_sign(value) == expected


def test_certified_sign_rejects_an_unseparated_ball() -> None:
    """A ball overlapping zero cannot be classified by its midpoint."""
    with pytest.raises(IntervalCertificationError):
        certified_sign(arb(0, 1))


def test_strictly_contains_distinguishes_boundary_contact() -> None:
    """Interior inclusion excludes a ball touching the outer boundary."""
    outer = arb(0, 2)
    assert strictly_contains(outer, arb(0, 1))
    assert not strictly_contains(arb(0, 1), arb(0, 1))


def test_interval_cholesky_certifies_spd_matrix() -> None:
    """The interval factor encloses a known positive-definite matrix."""
    matrix = arb_mat([[arb(4), arb(2)], [arb(2), arb(3)]])
    factor = interval_cholesky(matrix)
    assert (factor * factor.transpose()).contains(matrix)
    assert factor[0, 1].is_zero()
    assert factor[0, 0] > 0
    assert factor[1, 1] > 0


def test_interval_cholesky_rejects_uncertain_or_indefinite_pivot() -> None:
    """A pivot touching zero is not a positive-definiteness certificate."""
    with pytest.raises(IntervalCertificationError, match="pivot"):
        interval_cholesky(arb_mat([[arb(1), arb(0)], [arb(0), arb(0, 1)]]))
    with pytest.raises(IntervalCertificationError, match="pivot"):
        interval_cholesky(arb_mat([[arb(1), arb(2)], [arb(2), arb(1)]]))


def test_interval_cholesky_rejects_asymmetric_matrix() -> None:
    """The SPD contract does not silently symmetrize its input."""
    with pytest.raises(ValueError, match="symmetric"):
        interval_cholesky(arb_mat([[arb(1), arb(0)], [arb(1), arb(1)]]))


def test_krawczyk_image_certifies_linear_root() -> None:
    """An exact inverse contracts a linear two-variable system to its root."""
    jacobian = arb_mat([[arb(2), arb(0)], [arb(0), arb(4)]])
    preconditioner = arb_mat([[arb(1) / 2, arb(0)], [arb(0), arb(1) / 4]])
    center = (arb(1), arb(2))
    values = (arb(0), arb(0))
    domain = (arb(1, "0.25"), arb(2, "0.25"))
    image = krawczyk_image(center, values, jacobian, domain, preconditioner)
    box = arb_mat([[domain[0]], [domain[1]]])
    assert image[0, 0] == center[0]
    assert image[1, 0] == center[1]
    assert strictly_contains_vector(box, image)


def test_krawczyk_image_rejects_nonpoint_preconditioner() -> None:
    """An interval inverse cannot masquerade as the point preconditioner."""
    with pytest.raises(ValueError, match="preconditioner"):
        krawczyk_image(
            (arb(0),),
            (arb(0),),
            arb_mat([[arb(1)]]),
            (arb(0, 1),),
            arb_mat([[arb(1, "0.1")]]),
        )


def test_strictly_contains_vector_validates_dimensions() -> None:
    """Vector inclusion requires compatible interval column matrices."""
    with pytest.raises(ValueError, match="column"):
        strictly_contains_vector(arb_mat([[arb(0), arb(0)]]), arb_mat([[arb(0)]]))


def test_isolate_unique_root_proves_one_simple_root() -> None:
    """Interval Newton isolates the only root on the complete domain."""
    root = isolate_unique_root(lambda x: x * x - 2, lambda x: 2 * x, 1.0, 2.0)
    assert root.overlaps(arb(2).sqrt())
    assert float(root.upper() - root.lower()) < 1e-12
    assert (root * root - 2).contains(0)


def test_isolate_unique_root_rejects_multiple_roots() -> None:
    """A domain containing two roots is not silently reduced to one."""
    with pytest.raises(IntervalCertificationError, match="expected one"):
        isolate_unique_root(lambda x: x * x - 1, lambda x: 2 * x, -2.0, 2.0)


def test_isolate_unique_root_rejects_a_repeated_root() -> None:
    """A derivative interval containing zero cannot certify a simple root."""
    with pytest.raises(IntervalCertificationError):
        isolate_unique_root(lambda x: x * x, lambda x: 2 * x, -1.0, 1.0)


def test_best_root_float_selects_the_smaller_residual() -> None:
    """Adjacent binary64 candidates are ordered by certified residual bounds."""
    root = isolate_unique_root(lambda x: x * x - 2, lambda x: 2 * x, 1.0, 2.0)
    selected = best_root_float(lambda x: x * x - 2, root)
    neighbors = (
        math.nextafter(selected, -math.inf),
        math.nextafter(selected, math.inf),
    )
    assert selected == math.sqrt(2.0)
    assert abs(selected - math.sqrt(2.0)) <= min(
        abs(candidate - math.sqrt(2.0)) for candidate in neighbors
    )


def test_certified_float_accepts_a_unique_rounding_cell() -> None:
    """A narrow enclosure determines one native binary64 result."""
    value = math.pi
    with interval_precision(256):
        interval = ball_from_float(value) + arb(0, "1e-80")
        assert certified_float(interval) == value


def test_certified_float_rejects_multiple_binary64_candidates() -> None:
    """A wide enclosure must be refined rather than rounded speculatively."""
    with pytest.raises(IntervalCertificationError):
        certified_float(ball_from_bounds(1.0, math.nextafter(1.0, math.inf)))


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_nonfinite_binary64_input_is_rejected(bad: float) -> None:
    """Non-finite values never enter a geometric certificate."""
    with pytest.raises(ValueError):
        ball_from_float(bad)
