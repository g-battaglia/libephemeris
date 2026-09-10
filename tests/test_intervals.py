# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Tests for the certified Arb interval primitives."""

from __future__ import annotations

import math

import pytest
from flint import arb, ctx

from libephemeris.intervals import (
    IntervalCertificationError,
    ball_from_bounds,
    ball_from_float,
    certified_float,
    certified_sign,
    interval_precision,
    strictly_contains,
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
