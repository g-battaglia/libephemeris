# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Exact, disconnected checks for finite-source strict certificates.

This module checks caller-supplied rational geometry. It does not construct
certificates or obtain astronomical data. Its domain outcomes are deliberately
distinct from the physical reach and miss predicates.

Provenance:
    Project-authored exact-rational implementation of the fixed, reviewed
    finite-source strict-verifier specification in
    ``validation/golden/specs/occultation_contact_residuals/``
    (file ``finite-source-strict-verifier.md``; SHA-256
    e96d262c53b0b5af0edf70ca3f68f2a98f81987160ec4a8c84d69ee401faff46).
    The checker forms the specification's axis and cone residuals directly
    from supplied Fractions, classifies the ellipsoid by principal minors,
    and derives each dual affine form from evaluations at the origin and
    coordinate basis. No reference implementation, reference data, or
    astronomical values enter this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from fractions import Fraction
from typing import TypeAlias

_Vector: TypeAlias = tuple[Fraction, Fraction, Fraction]
_Matrix: TypeAlias = tuple[_Vector, _Vector, _Vector]
_ZERO = Fraction(0)
_ONE = Fraction(1)


class _DomainOutcome(Enum):
    """Reason a well-formed geometry is outside the proved checker domain."""

    UNSUPPORTED_SLOPE = auto()
    EXCLUDED_DOMAIN = auto()
    UNCERTIFIED_DOMAIN = auto()


class _OutsideDomain(ValueError):
    """Direct construction failed a well-formed theorem-domain check."""

    def __init__(self, outcome: _DomainOutcome) -> None:
        self.outcome = outcome
        super().__init__(outcome.name)


class _Target(Enum):
    """The three separate closed cone targets."""

    PENUMBRA = auto()
    UMBRA = auto()
    ANTUMBRA = auto()


@dataclass(frozen=True, slots=True)
class _DualPayload:
    """Exact coefficients of a target's dual-cone separator."""

    tau: Fraction
    beta: _Vector
    gate_multiplier: Fraction


def _scalar(value: object, name: str) -> Fraction:
    """Require an exact Fraction without coercing numerical substitutes."""
    if type(value) is not Fraction:
        raise TypeError(f"{name} must be an exact Fraction")
    return value


def _vector(value: object, name: str) -> _Vector:
    """Require exactly three Fraction coordinates in an immutable tuple."""
    if type(value) is not tuple or len(value) != 3:
        raise TypeError(f"{name} must be a three-tuple of Fractions")
    return (
        _scalar(value[0], f"{name}[0]"),
        _scalar(value[1], f"{name}[1]"),
        _scalar(value[2], f"{name}[2]"),
    )


def _matrix(value: object) -> _Matrix:
    """Require a three-row exact matrix before metric checks."""
    if type(value) is not tuple or len(value) != 3:
        raise TypeError("A must be a three-tuple of three-tuples")
    return (
        _vector(value[0], "A[0]"),
        _vector(value[1], "A[1]"),
        _vector(value[2], "A[2]"),
    )


def _dot(a: _Vector, b: _Vector) -> Fraction:
    """Compute a Euclidean dot product in exact rational arithmetic."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _subtract(a: _Vector, b: _Vector) -> _Vector:
    """Subtract two exact coordinate triples componentwise."""
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: _Vector, b: _Vector) -> _Vector:
    """Form the oriented cross product used by the cone residual."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _determinant(a: _Matrix) -> Fraction:
    """Expand the determinant of an exact three-by-three matrix."""
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def _principal_minors_nonnegative(a: _Matrix) -> bool:
    """Test all nonempty principal minors of a symmetric 3-by-3 matrix."""
    if any(a[i][i] < 0 for i in range(3)):
        return False
    for i, j in ((0, 1), (0, 2), (1, 2)):
        if a[i][i] * a[j][j] - a[i][j] * a[j][i] < 0:
            return False
    return _determinant(a) >= 0


def _inverse(a: _Matrix) -> _Matrix:
    """Invert an already certified positive-definite symmetric matrix."""
    det = _determinant(a)
    return (
        (
            (a[1][1] * a[2][2] - a[1][2] ** 2) / det,
            (a[0][2] * a[1][2] - a[0][1] * a[2][2]) / det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det,
        ),
        (
            (a[0][2] * a[1][2] - a[0][1] * a[2][2]) / det,
            (a[0][0] * a[2][2] - a[0][2] ** 2) / det,
            (a[0][1] * a[0][2] - a[0][0] * a[1][2]) / det,
        ),
        (
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det,
            (a[0][1] * a[0][2] - a[0][0] * a[1][2]) / det,
            (a[0][0] * a[1][1] - a[0][1] ** 2) / det,
        ),
    )


def _quadratic(a: _Matrix, x: _Vector) -> Fraction:
    """Evaluate the exact quadratic form x-transpose A x."""
    return sum((x[i] * _dot(a[i], x) for i in range(3)), _ZERO)


@dataclass(frozen=True, init=False, slots=True)
class _ValidatedGeometry:
    """Geometry whose entire solid receiver meets the reviewed domain.

    Even direct construction validates the domain and derives the axis and
    inverse. The checked values cannot be independently supplied.
    """

    B: _Vector
    p: _Vector
    r: Fraction
    m: Fraction
    A: _Matrix
    R: Fraction
    u: _Vector
    U: Fraction
    A_inverse: _Matrix

    def __init__(
        self, B: _Vector, p: _Vector, r: Fraction, m: Fraction, A: _Matrix, R: Fraction
    ) -> None:
        source = _vector(B, "B")
        occulter = _vector(p, "p")
        source_radius = _scalar(r, "r")
        occulter_radius = _scalar(m, "m")
        metric = _matrix(A)
        enclosure_radius = _scalar(R, "R")
        if source_radius < 0 or occulter_radius < 0:
            raise ValueError("r and m must be nonnegative")
        if enclosure_radius <= 0:
            raise ValueError("R must be positive")
        if any(metric[i][j] != metric[j][i] for i in range(3) for j in range(i + 1, 3)):
            raise ValueError("A must be symmetric")
        if (
            metric[0][0] <= 0
            or metric[0][0] * metric[1][1] - metric[0][1] ** 2 <= 0
            or _determinant(metric) <= 0
        ):
            raise ValueError("A must be positive definite")
        if occulter_radius == 0 or source_radius <= occulter_radius:
            raise _OutsideDomain(_DomainOutcome.UNSUPPORTED_SLOPE)
        axis = _subtract(occulter, source)
        axis_squared = _dot(axis, axis)
        if axis_squared <= (source_radius + occulter_radius) ** 2:
            raise _OutsideDomain(_DomainOutcome.EXCLUDED_DOMAIN)
        inverse_square = _ONE / (enclosure_radius * enclosure_radius)
        difference: _Matrix = tuple(
            tuple(
                metric[i][j] - (inverse_square if i == j else _ZERO) for j in range(3)
            )
            for i in range(3)
        )  # type: ignore[assignment]
        if not _principal_minors_nonnegative(difference):
            raise _OutsideDomain(_DomainOutcome.UNCERTIFIED_DOMAIN)
        if (
            _dot(source, source) <= (enclosure_radius + source_radius) ** 2
            or _dot(occulter, occulter) <= (enclosure_radius + occulter_radius) ** 2
        ):
            raise _OutsideDomain(_DomainOutcome.UNCERTIFIED_DOMAIN)
        object.__setattr__(self, "B", source)
        object.__setattr__(self, "p", occulter)
        object.__setattr__(self, "r", source_radius)
        object.__setattr__(self, "m", occulter_radius)
        object.__setattr__(self, "A", metric)
        object.__setattr__(self, "R", enclosure_radius)
        object.__setattr__(self, "u", axis)
        object.__setattr__(self, "U", axis_squared)
        object.__setattr__(self, "A_inverse", _inverse(metric))

    def _target_values(
        self, target: _Target, x: _Vector
    ) -> tuple[Fraction, Fraction | None, Fraction, _Vector]:
        v = _subtract(x, self.p)
        w = _dot(v, self.u)
        z = _cross(self.u, v)
        t = _dot(z, z)
        s = self.r + self.m
        d = self.r - self.m
        if target is _Target.PENUMBRA:
            return self.U * self.m + s * w, w + self.m * s, (self.U - s * s) * t, z
        if target is _Target.UMBRA:
            return self.U * self.m - d * w, w - self.m * d, (self.U - d * d) * t, z
        if target is _Target.ANTUMBRA:
            return d * w - self.U * self.m, None, (self.U - d * d) * t, z
        raise TypeError("target must be a private _Target")

    def contains(self, target: _Target, x: _Vector) -> bool:
        """Check closed target membership for an exact point in the solid K."""
        if type(target) is not _Target:
            raise TypeError("target must be a private _Target")
        point = _vector(x, "x")
        if _quadratic(self.A, point) > _ONE:
            raise ValueError("x lies outside K")
        a, g, h_t, _ = self._target_values(target, point)
        return (g is None or g >= 0) and a >= 0 and a * a >= h_t

    def strict_reach(self, target: _Target, x: _Vector) -> bool:
        """Validate an interior rational reach witness for one target."""
        if type(target) is not _Target:
            raise TypeError("target must be a private _Target")
        point = _vector(x, "x")
        if _quadratic(self.A, point) >= _ONE:
            return False
        a, g, h_t, _ = self._target_values(target, point)
        return (g is None or g > 0) and a > 0 and a * a > h_t

    def strict_miss(self, target: _Target, payload: _DualPayload) -> bool:
        """Validate a strict rational dual separator for one target."""
        if type(target) is not _Target:
            raise TypeError("target must be a private _Target")
        if type(payload) is not _DualPayload:
            raise TypeError("payload must be a private _DualPayload")
        tau = _scalar(payload.tau, "tau")
        beta = _vector(payload.beta, "beta")
        multiplier = _scalar(payload.gate_multiplier, "gate_multiplier")
        if tau < 0 or multiplier < 0:
            return False
        if target is _Target.ANTUMBRA and multiplier != 0:
            return False
        h = (
            self.U
            - (self.r + self.m if target is _Target.PENUMBRA else self.r - self.m) ** 2
        )
        if _dot(beta, beta) > tau * tau * h:
            return False

        def value(point: _Vector) -> Fraction:
            a, g, _, z = self._target_values(target, point)
            return (
                tau * a + _dot(beta, z) + (multiplier * g if g is not None else _ZERO)
            )

        origin = (_ZERO, _ZERO, _ZERO)
        c = value(origin)
        if c >= 0:
            return False
        n: _Vector = (
            value((_ONE, _ZERO, _ZERO)) - c,
            value((_ZERO, _ONE, _ZERO)) - c,
            value((_ZERO, _ZERO, _ONE)) - c,
        )
        return c * c > _quadratic(self.A_inverse, n)


def _prepare_geometry(
    B: _Vector, p: _Vector, r: Fraction, m: Fraction, A: _Matrix, R: Fraction
) -> _ValidatedGeometry | _DomainOutcome:
    """Validate exact inputs and classify the bounded theorem domain in order."""
    try:
        return _ValidatedGeometry(B, p, r, m, A, R)
    except _OutsideDomain as error:
        return error.outcome
