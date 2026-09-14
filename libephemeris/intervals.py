# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Certified arbitrary-precision interval primitives.

The geometry solvers use Arb real balls through :mod:`python-flint`. A ball is
an enclosure, not an estimate: every operation rounds outward and retains a
proof that the exact result lies inside. This module centralizes precision
management, exact binary64 ingestion, sign separation, and fail-closed
conversion back to a native float.

Arb is described by Fredrik Johansson, "Arb: efficient arbitrary-precision
midpoint-radius interval arithmetic", IEEE Transactions on Computers 66(8),
2017. The implementation is supplied by the MIT-licensed python-flint binding.

Provenance:
    Project-authored wrapper around the published Arb midpoint-radius interval
    arithmetic model and python-flint API. IEEE 754 defines the exact binary64
    bit decomposition used at the boundary. The module contains no astronomical
    coefficients and no values inferred from compatibility output.
"""

from __future__ import annotations

import math
import struct
from contextlib import contextmanager
from dataclasses import dataclass
from collections.abc import Callable, Sequence
from typing import Iterator, TypeAlias

from flint import arb, arb_mat, ctx

__all__ = [
    "Ball",
    "BallMatrix",
    "IntervalCertificationError",
    "ball_from_float",
    "ball_from_bounds",
    "best_root_float",
    "certified_float",
    "certified_sign",
    "contains_zero",
    "interval_cholesky",
    "interval_precision",
    "isolate_unique_root",
    "krawczyk_image",
    "strictly_contains",
    "strictly_contains_vector",
]

Ball: TypeAlias = arb
BallMatrix: TypeAlias = arb_mat


class IntervalCertificationError(ArithmeticError):
    """An interval calculation cannot yet prove the requested result."""


@dataclass(frozen=True, slots=True)
class _Binary64:
    """Exact rational decomposition of one finite binary64 value."""

    sign: int
    numerator: int
    denominator: int


def _binary64_ratio(value: float) -> _Binary64:
    """Return the exact signed rational represented by one finite binary64."""
    if not math.isfinite(value):
        raise ValueError("an interval input must be finite")
    bits = struct.unpack(">Q", struct.pack(">d", value))[0]
    sign = -1 if bits >> 63 else 1
    exponent = (bits >> 52) & 0x7FF
    fraction = bits & ((1 << 52) - 1)
    if exponent == 0:
        significand = fraction
        power = -1074
    else:
        significand = (1 << 52) | fraction
        power = int(exponent) - 1023 - 52
    if significand == 0:
        return _Binary64(sign, 0, 1)
    if power >= 0:
        return _Binary64(sign, significand << power, 1)
    return _Binary64(sign, significand, 1 << -power)


def ball_from_float(value: float) -> Ball:
    """Create an exact Arb ball from a finite Python float.

    Decimal formatting is not used: the binary64 numerator and denominator are
    reconstructed from their bits, so the ball contains exactly the value the
    caller supplied and no untracked input-rounding gap.
    """
    ratio = _binary64_ratio(value)
    result = arb(ratio.numerator) / arb(ratio.denominator)
    return -result if ratio.sign < 0 else result


def ball_from_bounds(lower: float, upper: float) -> Ball:
    """Enclose two finite binary64 bounds in one outward-rounded ball."""
    if lower > upper:
        raise ValueError("interval lower bound must not exceed upper bound")
    return ball_from_float(lower).union(ball_from_float(upper))


@contextmanager
def interval_precision(bits: int) -> Iterator[None]:
    """Temporarily set the process-local Arb working precision in bits."""
    if not isinstance(bits, int) or isinstance(bits, bool) or bits < 64:
        raise ValueError("interval precision must be an integer of at least 64 bits")
    previous = ctx.prec
    ctx.prec = bits
    try:
        yield
    finally:
        ctx.prec = previous


def contains_zero(value: Ball) -> bool:
    """Whether an Arb enclosure contains mathematical zero."""
    return bool(value.contains(0))


def certified_sign(value: Ball) -> int:
    """Return the proved sign of a ball, failing if it still contains zero."""
    if value > 0:
        return 1
    if value < 0:
        return -1
    if value.is_exact() and value.is_zero():
        return 0
    raise IntervalCertificationError("interval sign is not separated from zero")


def strictly_contains(outer: Ball, inner: Ball) -> bool:
    """Whether ``inner`` is proved to lie in the interior of ``outer``."""
    return bool(outer.contains_interior(inner))


def strictly_contains_vector(outer: BallMatrix, inner: BallMatrix) -> bool:
    """Whether one interval column vector strictly contains another.

    Args:
        outer: An ``n x 1`` interval vector.
        inner: An ``n x 1`` interval vector of the same dimension.

    Returns:
        True only when every inner component lies in the interior of its
        corresponding outer component.

    Raises:
        ValueError: If either matrix is not a compatible column vector.
    """
    if outer.ncols() != 1 or inner.ncols() != 1 or outer.nrows() != inner.nrows():
        raise ValueError("interval vectors must be compatible column matrices")
    return all(
        strictly_contains(outer[index, 0], inner[index, 0])
        for index in range(outer.nrows())
    )


def interval_cholesky(matrix: BallMatrix) -> BallMatrix:
    """Certify positive definiteness and return a lower Cholesky enclosure.

    The standard scalar Cholesky recurrence is evaluated in Arb. Every pivot
    must be proved strictly positive before its square root is taken; an
    enclosure touching zero is not accepted as numerical positive definiteness.

    Args:
        matrix: Symmetric square interval matrix.

    Returns:
        Lower-triangular interval matrix ``L`` enclosing a factor with
        ``matrix = L * L.T``.

    Raises:
        ValueError: If the matrix is empty, non-square, or not certifiably
            symmetric.
        IntervalCertificationError: If a pivot is not proved positive.
    """
    rows = matrix.nrows()
    if rows == 0 or matrix.ncols() != rows:
        raise ValueError("Cholesky input must be a non-empty square matrix")
    for row in range(rows):
        for column in range(row):
            if matrix[row, column] != matrix[column, row]:
                raise ValueError("Cholesky input must be certifiably symmetric")
    factor = arb_mat(rows, rows)
    for row in range(rows):
        for column in range(row + 1):
            remainder = matrix[row, column]
            for index in range(column):
                remainder -= factor[row, index] * factor[column, index]
            if row == column:
                if not remainder > 0:
                    raise IntervalCertificationError(
                        f"Cholesky pivot {row} is not separated above zero"
                    )
                factor[row, column] = remainder.sqrt()
            else:
                factor[row, column] = remainder / factor[column, column]
    return factor


def _column_vector(values: Sequence[Ball], name: str) -> BallMatrix:
    """Build a column matrix after validating one interval-vector input."""
    if not values:
        raise ValueError(f"{name} must be a non-empty interval vector")
    return arb_mat([[value] for value in values])


def krawczyk_image(
    center: Sequence[Ball],
    values_at_center: Sequence[Ball],
    jacobian_box: BallMatrix,
    domain: Sequence[Ball],
    preconditioner: BallMatrix,
) -> BallMatrix:
    """Return the interval Krawczyk image for a square nonlinear system.

    For a point ``x`` and interval box ``X``, system value ``F(x)``, interval
    Jacobian ``J(X)``, and point preconditioner ``C``, the operator is

    ``K = x - C F(x) + (I - C J(X)) (X - x)``.

    Strict inclusion ``K`` inside ``X`` proves one root in the box under the
    usual Krawczyk theorem. This helper computes only the outward enclosure;
    callers retain responsibility for complete-domain subdivision and every
    problem-specific boundary or active-set certificate.

    Args:
        center: Exact or point-valued center vector ``x``.
        values_at_center: Outward enclosure of ``F(x)``.
        jacobian_box: Outward enclosure of ``J(X)``.
        domain: Interval box ``X``.
        preconditioner: Square point matrix ``C``.

    Returns:
        The outward-rounded Krawczyk image as an ``n x 1`` matrix.

    Raises:
        ValueError: If dimensions do not form one square system, inputs are not
            finite, point quantities are not exact, the preconditioner is
            singular, or the center does not lie inside the domain.
    """
    dimension = len(center)
    if dimension == 0 or len(values_at_center) != dimension or len(domain) != dimension:
        raise ValueError("Krawczyk vectors must have the same non-zero dimension")
    if jacobian_box.nrows() != dimension or jacobian_box.ncols() != dimension:
        raise ValueError("Krawczyk Jacobian must be square and match the vectors")
    if preconditioner.nrows() != dimension or preconditioner.ncols() != dimension:
        raise ValueError("Krawczyk preconditioner must match the system dimension")
    if any(not value.is_exact() or not value.is_finite() for value in center):
        raise ValueError("Krawczyk center must contain finite exact point values")
    if any(not value.is_finite() for value in values_at_center):
        raise ValueError("Krawczyk function values must be finite intervals")
    if any(not value.is_finite() for value in domain):
        raise ValueError("Krawczyk domain must contain finite intervals")
    if any(
        not jacobian_box[row, column].is_finite()
        for row in range(dimension)
        for column in range(dimension)
    ):
        raise ValueError("Krawczyk Jacobian must contain finite intervals")
    if any(
        not preconditioner[row, column].is_exact()
        or not preconditioner[row, column].is_finite()
        for row in range(dimension)
        for column in range(dimension)
    ):
        raise ValueError("Krawczyk preconditioner must be a finite exact point matrix")
    if not preconditioner.det().is_finite() or contains_zero(preconditioner.det()):
        raise ValueError("Krawczyk preconditioner must be nonsingular")
    if any(not domain[index].contains(center[index]) for index in range(dimension)):
        raise ValueError("Krawczyk center must lie inside the domain box")

    center_vector = _column_vector(center, "center")
    value_vector = _column_vector(values_at_center, "values_at_center")
    domain_vector = _column_vector(domain, "domain")
    identity = arb_mat(
        [
            [arb(1 if row == column else 0) for column in range(dimension)]
            for row in range(dimension)
        ]
    )
    return (
        center_vector
        - preconditioner * value_vector
        + (identity - preconditioner * jacobian_box) * (domain_vector - center_vector)
    )


def _exact_lower(value: Ball) -> Ball:
    """Return the exact dyadic lower endpoint of an Arb ball."""
    return arb(value.lower().mid())


def _exact_upper(value: Ball) -> Ball:
    """Return the exact dyadic upper endpoint of an Arb ball."""
    return arb(value.upper().mid())


def isolate_unique_root(
    function: Callable[[Ball], Ball],
    derivative: Callable[[Ball], Ball],
    lower: float,
    upper: float,
    *,
    precision: int = 192,
    max_boxes: int = 1_000_000,
) -> Ball:
    """Isolate and certify one simple root over a complete closed interval.

    Every discarded box has a function enclosure excluding zero. A retained
    box is contracted with interval Newton; uniqueness is proved only when its
    derivative enclosure excludes zero and the Newton image lies strictly in
    the box. The entire input interval is exhausted, so a second root cannot be
    hidden outside the returned enclosure.

    Args:
        function: Outward-rounded interval extension of a continuous scalar.
        derivative: Outward-rounded interval extension of its derivative.
        lower: Finite lower endpoint.
        upper: Finite upper endpoint, strictly greater than ``lower``.
        precision: Arb precision in bits.
        max_boxes: Hard arithmetic work limit; exceeding it fails closed.

    Returns:
        An Arb enclosure containing exactly one simple root.

    Raises:
        ValueError: If the interval or limits are invalid.
        IntervalCertificationError: If existence or uniqueness is not proved.
    """
    if not math.isfinite(lower) or not math.isfinite(upper) or not lower < upper:
        raise ValueError("root interval must have finite increasing endpoints")
    if not isinstance(max_boxes, int) or isinstance(max_boxes, bool) or max_boxes <= 0:
        raise ValueError("max_boxes must be a positive integer")

    with interval_precision(precision):
        left = ball_from_float(lower)
        right = ball_from_float(upper)
        queue: list[tuple[Ball, Ball]] = [(left, right)]
        roots: list[Ball] = []
        boxes = 0
        while queue:
            lo, hi = queue.pop()
            boxes += 1
            if boxes > max_boxes:
                raise IntervalCertificationError(
                    "interval subdivision exceeded its certified work limit"
                )
            domain = lo.union(hi)
            values = function(domain)
            if not contains_zero(values):
                continue

            slope = derivative(domain)
            midpoint = (lo + hi) / 2
            point_value = function(midpoint)
            # An exact interior zero is already a certified singleton root when
            # the derivative excludes zero.  Handling it before interval Newton
            # avoids repeatedly subdividing an exact dyadic solution.
            if (
                point_value.is_exact()
                and point_value.is_zero()
                and not contains_zero(slope)
            ):
                roots.append(midpoint)
                continue
            if not slope.is_finite():
                if midpoint == lo or midpoint == hi:
                    raise IntervalCertificationError(
                        "root derivative interval is not finite"
                    )
                queue.append((lo, midpoint))
                queue.append((midpoint, hi))
                continue
            if not contains_zero(slope):
                newton = midpoint - point_value / slope
                try:
                    contracted = domain.intersection(newton)
                except ValueError:
                    continue
                if contracted.is_finite():
                    contracted_lo = _exact_lower(contracted)
                    contracted_hi = _exact_upper(contracted)
                    if contracted_lo < contracted_hi:
                        if strictly_contains(domain, contracted):
                            domain = contracted
                            lo, hi = contracted_lo, contracted_hi
                        left_sign = certified_sign(function(lo))
                        right_sign = certified_sign(function(hi))
                        if left_sign * right_sign < 0 and strictly_contains(
                            lo.union(hi), midpoint - point_value / derivative(domain)
                        ):
                            certified = lo.union(hi)
                            for _ in range(32):
                                center = certified.mid()
                                try:
                                    narrowed = certified.intersection(
                                        center
                                        - function(center) / derivative(certified)
                                    )
                                except ValueError:
                                    break
                                if narrowed == certified:
                                    break
                                certified = narrowed
                            roots.append(certified)
                            continue

            midpoint = (lo + hi) / 2
            if midpoint == lo or midpoint == hi:
                raise IntervalCertificationError(
                    "root box reached the arithmetic subdivision floor"
                )
            queue.append((lo, midpoint))
            queue.append((midpoint, hi))

        if len(roots) != 1:
            raise IntervalCertificationError(
                f"expected one certified root, isolated {len(roots)}"
            )
        return roots[0]


def best_root_float(function: Callable[[Ball], Ball], root: Ball) -> float:
    """Select the binary64 candidate with the smaller certified residual.

    The two adjacent candidates around the isolated root are evaluated as exact
    binary64 inputs. One is selected only when its absolute residual enclosure
    is strictly below the other's. An exact tie uses the earlier number.
    """
    midpoint = float(root.mid())
    candidates = tuple(
        sorted(
            {
                math.nextafter(midpoint, -math.inf),
                midpoint,
                math.nextafter(midpoint, math.inf),
            }
        )
    )
    for bits in (192, 256, 384, 512):
        with interval_precision(bits):
            residuals = [
                (candidate, abs(function(ball_from_float(candidate))))
                for candidate in candidates
            ]
            ordered = sorted(residuals, key=lambda item: float(item[1].mid()))
            best_candidate, best_residual = ordered[0]
            second_residual = ordered[1][1]
            if best_residual.upper() < second_residual.lower():
                return best_candidate
            if (
                best_residual.is_exact()
                and second_residual.is_exact()
                and best_residual == second_residual
            ):
                return min(
                    candidate
                    for candidate, residual in residuals
                    if residual == best_residual
                )
    raise IntervalCertificationError(
        "adjacent binary64 root candidates have overlapping residual enclosures"
    )


def certified_float(value: Ball) -> float:
    """Return the unique nearest binary64 enclosed by a sufficiently narrow ball.

    The midpoint is rounded by Python to binary64, then the two adjacent
    representable values define its rounding cell. The result is accepted only
    when the complete Arb enclosure lies strictly inside that cell. Otherwise
    callers must increase precision rather than guess one side.
    """
    if not value.is_finite():
        raise IntervalCertificationError("interval result is not finite")
    candidate = float(value.mid())
    if not math.isfinite(candidate):
        raise IntervalCertificationError("interval midpoint is outside binary64")
    lower_neighbor = math.nextafter(candidate, -math.inf)
    upper_neighbor = math.nextafter(candidate, math.inf)
    lower_boundary = (ball_from_float(lower_neighbor) + ball_from_float(candidate)) / 2
    upper_boundary = (ball_from_float(candidate) + ball_from_float(upper_neighbor)) / 2
    if value.lower() > lower_boundary and value.upper() < upper_boundary:
        return candidate
    if value.is_exact() and value == ball_from_float(candidate):
        return candidate
    raise IntervalCertificationError(
        "interval does not determine a unique nearest binary64 value"
    )
