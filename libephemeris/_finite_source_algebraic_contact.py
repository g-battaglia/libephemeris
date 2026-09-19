# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disconnected exact checks for finite-source algebraic contact certificates.

The caller supplies a rational geometry and one algebraic dual payload. This
module neither constructs a certificate nor reads astronomical state. It is
intentionally unbounded and must not be connected to a runtime producer
without a separately reviewed limit covering the entire invocation.

Provenance:
    Project-authored implementation of the private specification
    ``validation/golden/specs/occultation_contact_residuals/``
    ``finite-source-algebraic-contact-verifier.md`` (SHA-256
    b2aca7a287a9fb41711e2b6fe74c2cdc1f92bc396c4fa05934e722bb4a5508cd).
    All geometry, polynomial and sign operations are derived from that
    specification and caller-supplied exact values.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import TypeAlias

from libephemeris._finite_source_certificates import _Target, _ValidatedGeometry

_Poly: TypeAlias = tuple[Fraction, ...]
_Vector: TypeAlias = tuple[Fraction, Fraction, Fraction]
_Matrix: TypeAlias = tuple[_Vector, _Vector, _Vector]
_ZERO = Fraction(0)
_ONE = Fraction(1)
_POLY_ZERO: _Poly = (_ZERO,)


@dataclass(frozen=True, slots=True)
class _AlgebraicContactPayload:
    """One isolated real root and dual coefficients in its rational field."""

    P: tuple[int, ...]
    lo: Fraction
    hi: Fraction
    tau: _Poly
    beta: tuple[_Poly, _Poly, _Poly]
    gate_multiplier: _Poly | None


@dataclass(frozen=True, slots=True)
class _AffineTarget:
    """Exact affine cone and optional physical gate for one target."""

    a0: Fraction
    v: _Vector
    z0: _Vector
    Z: _Matrix
    g0: Fraction | None
    t: _Vector | None
    H: Fraction


def _trim(coefficients: list[Fraction]) -> _Poly:
    """Return a canonical ascending-coefficient rational polynomial."""
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients) if coefficients else _POLY_ZERO


def _add(left: _Poly, right: _Poly) -> _Poly:
    """Add exact rational polynomials."""
    return _trim(
        [
            (left[i] if i < len(left) else _ZERO)
            + (right[i] if i < len(right) else _ZERO)
            for i in range(max(len(left), len(right)))
        ]
    )


def _scale(poly: _Poly, factor: Fraction) -> _Poly:
    """Multiply a polynomial by an exact rational scalar."""
    return _trim([coefficient * factor for coefficient in poly])


def _sub(left: _Poly, right: _Poly) -> _Poly:
    """Subtract exact rational polynomials."""
    return _add(left, _scale(right, -_ONE))


def _mul(left: _Poly, right: _Poly) -> _Poly:
    """Multiply exact rational polynomials."""
    if left == _POLY_ZERO or right == _POLY_ZERO:
        return _POLY_ZERO
    result = [_ZERO] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return _trim(result)


def _derivative(poly: _Poly) -> _Poly:
    """Differentiate an exact rational polynomial."""
    return _trim([poly[i] * i for i in range(1, len(poly))])


def _remainder(dividend: _Poly, divisor: _Poly) -> _Poly:
    """Return the exact Euclidean polynomial remainder."""
    if divisor == _POLY_ZERO:
        raise ZeroDivisionError("zero polynomial divisor")
    rest = dividend
    while rest != _POLY_ZERO and len(rest) >= len(divisor):
        factor = rest[-1] / divisor[-1]
        shift = len(rest) - len(divisor)
        subtraction = [_ZERO] * shift + [factor * x for x in divisor]
        rest = _sub(rest, tuple(subtraction))
    return rest


def _gcd(left: _Poly, right: _Poly) -> _Poly:
    """Compute the monic rational-polynomial greatest common divisor."""
    while right != _POLY_ZERO:
        left, right = right, _remainder(left, right)
    if left == _POLY_ZERO:
        return _POLY_ZERO
    return _scale(left, _ONE / left[-1])


def _evaluate(poly: _Poly, point: Fraction) -> Fraction:
    """Evaluate an exact rational polynomial by Horner's rule."""
    result = _ZERO
    for coefficient in reversed(poly):
        result = result * point + coefficient
    return result


def _sturm_chain(poly: _Poly) -> tuple[_Poly, ...]:
    """Build the signed Sturm chain of a squarefree polynomial."""
    if len(poly) == 1:
        return (poly,)
    chain = [poly, _derivative(poly)]
    while chain[-1] != _POLY_ZERO and len(chain[-1]) > 1:
        rest = _remainder(chain[-2], chain[-1])
        if rest == _POLY_ZERO:
            raise ValueError("Sturm input must be squarefree")
        chain.append(_scale(rest, -_ONE))
    return tuple(chain)


def _variation(chain: tuple[_Poly, ...], point: Fraction) -> int:
    """Count sign changes after dropping zero endpoint evaluations."""
    signs = []
    for poly in chain:
        value = _evaluate(poly, point)
        if value != 0:
            signs.append(1 if value > 0 else -1)
    return sum(signs[i] != signs[i + 1] for i in range(len(signs) - 1))


def _root_count(chain: tuple[_Poly, ...], lo: Fraction, hi: Fraction) -> int:
    """Count distinct roots in an open interval with nonroot endpoints."""
    return _variation(chain, lo) - _variation(chain, hi)


def _interval_multiply(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    """Multiply two exact closed rational intervals."""
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def _interval_evaluate(
    poly: _Poly, lo: Fraction, hi: Fraction
) -> tuple[Fraction, Fraction]:
    """Enclose all polynomial values on a rational interval by Horner."""
    result = (_ZERO, _ZERO)
    domain = (lo, hi)
    for coefficient in reversed(poly):
        low, high = _interval_multiply(result, domain)
        result = (low + coefficient, high + coefficient)
    return result


class _RootSign:
    """Decide exact signs at one isolated real algebraic root."""

    def __init__(self, poly: _Poly, lo: Fraction, hi: Fraction) -> None:
        self.poly = poly
        self.lo = lo
        self.hi = hi
        self.sturm = _sturm_chain(poly)

    def sign(self, expression: _Poly) -> int:
        """Return -1, 0 or 1 for an exact polynomial at the selected root."""
        if expression == _POLY_ZERO:
            return 0
        common = _gcd(self.poly, expression)
        if len(common) > 1 and _root_count(_sturm_chain(common), self.lo, self.hi) == 1:
            return 0

        lo, hi = self.lo, self.hi
        while True:
            lower, upper = _interval_evaluate(expression, lo, hi)
            if lower > 0:
                return 1
            if upper < 0:
                return -1
            midpoint = (lo + hi) / 2
            if _evaluate(self.poly, midpoint) == 0:
                value = _evaluate(expression, midpoint)
                return (value > 0) - (value < 0)
            if _root_count(self.sturm, lo, midpoint) == 1:
                hi = midpoint
            else:
                lo = midpoint


def _check_value_poly(value: object, name: str) -> _Poly:
    """Require a nonempty, trimmed tuple of exact Fraction coefficients."""
    if type(value) is not tuple or not value:
        raise TypeError(f"{name} must be a nonempty polynomial tuple")
    for coefficient in value:
        if type(coefficient) is not Fraction:
            raise TypeError(f"{name} coefficients must be exact Fractions")
    if len(value) > 1 and value[-1] == 0:
        raise ValueError(f"{name} must not have a trailing zero")
    return value


def _check_structure(
    geometry: object, target: object, payload: object
) -> tuple[_ValidatedGeometry, _Target, _AlgebraicContactPayload]:
    """Check only object, tuple and coefficient syntax before domain work."""
    if type(geometry) is not _ValidatedGeometry:
        raise TypeError("geometry must be an exact _ValidatedGeometry")
    if type(target) is not _Target:
        raise TypeError("target must be an exact _Target")
    if type(payload) is not _AlgebraicContactPayload:
        raise TypeError("payload must be an exact _AlgebraicContactPayload")
    if type(payload.P) is not tuple or not payload.P:
        raise TypeError("P must be a nonempty tuple of integers")
    if any(type(coefficient) is not int for coefficient in payload.P):
        raise TypeError("P coefficients must be exact integers")
    if type(payload.lo) is not Fraction or type(payload.hi) is not Fraction:
        raise TypeError("root interval endpoints must be exact Fractions")
    _check_value_poly(payload.tau, "tau")
    if type(payload.beta) is not tuple or len(payload.beta) != 3:
        raise TypeError("beta must be a three-tuple of polynomials")
    for index, beta in enumerate(payload.beta):
        _check_value_poly(beta, f"beta[{index}]")
    if target is _Target.ANTUMBRA:
        if payload.gate_multiplier is not None:
            raise TypeError("antumbra has no gate multiplier")
    else:
        _check_value_poly(payload.gate_multiplier, "gate_multiplier")
    return geometry, target, payload


def _check_root(payload: _AlgebraicContactPayload) -> _RootSign | None:
    """Return a sign engine only for a normalized one-root certificate."""
    coefficients = payload.P
    if len(coefficients) < 2 or coefficients[-1] <= 0:
        return None
    content = 0
    for coefficient in coefficients:
        content = gcd(content, abs(coefficient))
    if content != 1 or payload.lo >= payload.hi:
        return None
    poly = tuple(Fraction(coefficient) for coefficient in coefficients)
    if len(_gcd(poly, _derivative(poly))) > 1:
        return None
    if _evaluate(poly, payload.lo) == 0 or _evaluate(poly, payload.hi) == 0:
        return None
    engine = _RootSign(poly, payload.lo, payload.hi)
    return engine if _root_count(engine.sturm, payload.lo, payload.hi) == 1 else None


def _dot(left: _Vector, right: _Vector) -> Fraction:
    """Dot two rational three-vectors."""
    return sum((left[i] * right[i] for i in range(3)), _ZERO)


def _cross(left: _Vector, right: _Vector) -> _Vector:
    """Form the oriented cross product of two rational three-vectors."""
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _vector_scale(vector: _Vector, factor: Fraction) -> _Vector:
    """Scale exactly three rational coordinates."""
    return (factor * vector[0], factor * vector[1], factor * vector[2])


def _affine_target(geometry: _ValidatedGeometry, target: _Target) -> _AffineTarget:
    """Derive the physical target's exact affine forms from source geometry."""
    u = geometry.u
    U = geometry.U
    s = geometry.r + geometry.m
    d = geometry.r - geometry.m
    w0 = -_dot(u, geometry.p)
    z0 = _vector_scale(_cross(u, geometry.p), -_ONE)
    Z: _Matrix = (
        (_ZERO, -u[2], u[1]),
        (u[2], _ZERO, -u[0]),
        (-u[1], u[0], _ZERO),
    )
    if target is _Target.PENUMBRA:
        return _AffineTarget(
            U * geometry.m + s * w0,
            _vector_scale(u, s),
            z0,
            Z,
            w0 + geometry.m * s,
            u,
            U - s * s,
        )
    if target is _Target.UMBRA:
        return _AffineTarget(
            U * geometry.m - d * w0,
            _vector_scale(u, -d),
            z0,
            Z,
            w0 - geometry.m * d,
            u,
            U - d * d,
        )
    return _AffineTarget(
        d * w0 - U * geometry.m,
        _vector_scale(u, d),
        z0,
        Z,
        None,
        None,
        U - d * d,
    )


def _sum_scaled(polys: tuple[_Poly, ...], factors: tuple[Fraction, ...]) -> _Poly:
    """Form one exact rational linear combination of polynomials."""
    result = _POLY_ZERO
    for poly, factor in zip(polys, factors, strict=True):
        result = _add(result, _scale(poly, factor))
    return result


def _gate_parts(affine: _AffineTarget) -> tuple[Fraction, _Vector]:
    """Require the internally derived gate for a gated target."""
    if affine.g0 is None or affine.t is None:
        raise RuntimeError("gated target is missing its derived gate")
    return affine.g0, affine.t


def _verify_algebraic_contact(
    geometry: _ValidatedGeometry,
    target: _Target,
    payload: _AlgebraicContactPayload,
) -> bool:
    """Check one exact algebraic dual certificate for per-target contact.

    Args:
        geometry: Caller-supplied private geometry, revalidated from base fields.
        target: One physical cone target.
        payload: One isolated root and dual polynomials in that root.

    Returns:
        Whether every exact contact-certificate condition passes.

    Raises:
        TypeError: An object or coefficient has the wrong exact type.
        ValueError: A structural polynomial encoding or geometry is malformed.
        _OutsideDomain: The reconstructed geometry fails the theorem domain.
    """
    geometry, target, payload = _check_structure(geometry, target, payload)
    fresh = _ValidatedGeometry(
        geometry.B, geometry.p, geometry.r, geometry.m, geometry.A, geometry.R
    )
    sign_engine = _check_root(payload)
    if sign_engine is None:
        return False
    degree = len(payload.P) - 1
    values: list[_Poly] = [payload.tau, *payload.beta]
    if target is not _Target.ANTUMBRA:
        if payload.gate_multiplier is None:
            raise TypeError("gated target requires a gate multiplier")
        values.append(payload.gate_multiplier)
    if any(len(poly) > degree for poly in values):
        raise ValueError("value polynomial degree must be below root degree")

    affine = _affine_target(fresh, target)
    tau = payload.tau
    beta = payload.beta
    gamma = payload.gate_multiplier
    if sign_engine.sign(tau) < 0:
        return False
    if gamma is not None and sign_engine.sign(gamma) < 0:
        return False
    dual_slack = _sub(
        _scale(_mul(tau, tau), affine.H),
        _add(
            _add(_mul(beta[0], beta[0]), _mul(beta[1], beta[1])), _mul(beta[2], beta[2])
        ),
    )
    if sign_engine.sign(dual_slack) < 0:
        return False

    normal = []
    for i in range(3):
        component = _add(
            _scale(tau, affine.v[i]),
            _sum_scaled(beta, tuple(affine.Z[j][i] for j in range(3))),
        )
        if gamma is not None:
            _, gate_t = _gate_parts(affine)
            component = _add(component, _scale(gamma, gate_t[i]))
        normal.append(component)
    n = tuple(normal)
    c = _add(_scale(tau, affine.a0), _sum_scaled(beta, affine.z0))
    if gamma is not None:
        gate_g0, _ = _gate_parts(affine)
        c = _add(c, _scale(gamma, gate_g0))
    if sign_engine.sign(c) >= 0:
        return False

    b = tuple(_sum_scaled(n, fresh.A_inverse[i]) for i in range(3))
    support = _sub(
        _mul(c, c),
        _add(_add(_mul(n[0], b[0]), _mul(n[1], b[1])), _mul(n[2], b[2])),
    )
    if sign_engine.sign(support) != 0:
        return False

    a_num = _sub(_scale(c, affine.a0), _sum_scaled(b, affine.v))
    if sign_engine.sign(a_num) > 0:
        return False
    z_num = tuple(
        _sub(_scale(c, affine.z0[j]), _sum_scaled(b, affine.Z[j])) for j in range(3)
    )
    cone_slack = _sub(
        _mul(a_num, a_num),
        _scale(
            _add(
                _add(_mul(z_num[0], z_num[0]), _mul(z_num[1], z_num[1])),
                _mul(z_num[2], z_num[2]),
            ),
            affine.H,
        ),
    )
    if sign_engine.sign(cone_slack) < 0:
        return False
    if affine.g0 is not None:
        gate_g0, gate_t = _gate_parts(affine)
        g_num = _sub(_scale(c, gate_g0), _sum_scaled(b, gate_t))
        if sign_engine.sign(g_num) > 0:
            return False
    return True
