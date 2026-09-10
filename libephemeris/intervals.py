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
from typing import Iterator, TypeAlias

from flint import arb, arb_mat, ctx

__all__ = [
    "Ball",
    "BallMatrix",
    "IntervalCertificationError",
    "ball_from_float",
    "ball_from_bounds",
    "certified_float",
    "certified_sign",
    "contains_zero",
    "interval_precision",
    "strictly_contains",
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
