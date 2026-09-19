# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Exact real-algebraic arithmetic for the disconnected contact generator.

This producer-side module uses rational FLINT polynomials for factorization
and field operations. It has no dependency on the contact verifier's sign
engine. Every real-root choice and comparison is certified with rational
Sturm counts and exact interval arithmetic; no ball midpoint is decisive.

Provenance:
    Project-authored implementation of the private generator specification
    ``validation/golden/specs/occultation_contact_residuals/``
    ``finite-source-contact-generator.md`` (SHA-256
    9ea40fb03e873f4149985a09010cebdd5c8427d759c830bd4998b37fa9270b7c).
    Exact rational factorization and matrices use the project's existing
    python-flint dependency; the contact verifier remains independent.
"""

from __future__ import annotations

from fractions import Fraction
from math import floor, gcd, lcm
from typing import TypeAlias

from flint import fmpq, fmpq_poly

_Poly: TypeAlias = tuple[Fraction, ...]
_Interval: TypeAlias = tuple[Fraction, Fraction]
_ZERO = Fraction(0)
_ONE = Fraction(1)


def _as_flint(value: Fraction | int) -> fmpq:
    """Convert one exact Python rational without a decimal intermediary."""
    rational = Fraction(value)
    return fmpq(rational.numerator, rational.denominator)


def _as_fraction(value: fmpq) -> Fraction:
    """Convert one FLINT rational to a normalized Python Fraction."""
    return Fraction(int(value.p), int(value.q))


def _poly_coefficients(poly: fmpq_poly) -> _Poly:
    """Return ascending exact rational coefficients, including one zero."""
    coefficients = tuple(_as_fraction(item) for item in poly.coeffs())
    return coefficients or (_ZERO,)


def _integer_polynomial(poly: fmpq_poly) -> tuple[int, ...]:
    """Normalize a nonconstant rational polynomial to primitive integers."""
    coefficients = _poly_coefficients(poly)
    if len(coefficients) < 2:
        raise ValueError("a root polynomial must be nonconstant")
    denominator = lcm(*(item.denominator for item in coefficients))
    integers = tuple(int(item * denominator) for item in coefficients)
    content = gcd(*(abs(item) for item in integers))
    if content == 0:
        raise ValueError("zero root polynomial")
    sign = 1 if integers[-1] > 0 else -1
    return tuple(sign * item // content for item in integers)


def _flint_polynomial(coefficients: _Poly | tuple[int, ...]) -> fmpq_poly:
    """Create an exact FLINT polynomial from ascending coefficients."""
    return fmpq_poly([_as_flint(item) for item in coefficients])


def _interval_add(left: _Interval, right: _Interval) -> _Interval:
    """Add closed rational intervals exactly."""
    return left[0] + right[0], left[1] + right[1]


def _interval_neg(value: _Interval) -> _Interval:
    """Negate a closed rational interval."""
    return -value[1], -value[0]


def _interval_sub(left: _Interval, right: _Interval) -> _Interval:
    """Subtract closed rational intervals exactly."""
    return _interval_add(left, _interval_neg(right))


def _interval_mul(left: _Interval, right: _Interval) -> _Interval:
    """Multiply closed rational intervals by checking all four corners."""
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def _interval_div(left: _Interval, right: _Interval) -> _Interval:
    """Divide by an interval certified away from zero."""
    if right[0] <= 0 <= right[1]:
        raise ZeroDivisionError("interval divisor contains zero")
    reciprocal = (1 / right[1], 1 / right[0])
    return _interval_mul(left, reciprocal)


def _interval_evaluate(poly: fmpq_poly, interval: _Interval) -> _Interval:
    """Enclose all values of one rational polynomial by interval Horner."""
    result = (_ZERO, _ZERO)
    for coefficient in reversed(_poly_coefficients(poly)):
        result = _interval_add(_interval_mul(result, interval), (coefficient,) * 2)
    return result


def _sturm_chain(poly: fmpq_poly) -> tuple[fmpq_poly, ...]:
    """Build an exact signed Sturm chain for a squarefree polynomial."""
    if poly.degree() < 1:
        raise ValueError("Sturm input must be nonconstant")
    chain = [poly, poly.derivative()]
    while chain[-1].degree() > 0:
        remainder = chain[-2] % chain[-1]
        if remainder == 0:
            raise ValueError("Sturm input must be squarefree")
        chain.append(-remainder)
    return tuple(chain)


def _variation(chain: tuple[fmpq_poly, ...], point: Fraction) -> int:
    """Count Sturm sign changes after deleting zero endpoint values."""
    argument = _as_flint(point)
    signs: list[int] = []
    for polynomial in chain:
        value = polynomial(argument)
        if value != 0:
            signs.append(1 if value > 0 else -1)
    return sum(signs[i] != signs[i + 1] for i in range(len(signs) - 1))


def _root_count(chain: tuple[fmpq_poly, ...], lo: Fraction, hi: Fraction) -> int:
    """Count roots in an open interval with nonroot endpoints."""
    if lo >= hi:
        raise ValueError("root interval must increase")
    if chain[0](_as_flint(lo)) == 0 or chain[0](_as_flint(hi)) == 0:
        raise ValueError("root interval endpoint is a root")
    return _variation(chain, lo) - _variation(chain, hi)


class _RealRoot:
    """One exact real root of an irreducible rational polynomial."""

    __slots__ = ("P", "poly", "sturm", "lo", "hi", "_canonical")

    def __init__(self, P: tuple[int, ...], lo: Fraction, hi: Fraction) -> None:
        self.P = P
        self.poly = _flint_polynomial(P)
        self.sturm = _sturm_chain(self.poly)
        if _root_count(self.sturm, lo, hi) != 1:
            raise ValueError("interval must isolate one root")
        self.lo = lo
        self.hi = hi
        self._canonical: _Interval | None = None

    def refine(self) -> None:
        """Bisect the isolator without changing its selected root."""
        midpoint = (self.lo + self.hi) / 2
        if self.poly(_as_flint(midpoint)) == 0:
            # Only a degree-one minimal polynomial can have a rational root.
            width = self.hi - self.lo
            self.lo = midpoint - width / 4
            self.hi = midpoint + width / 4
            return
        if _root_count(self.sturm, self.lo, midpoint) == 1:
            self.hi = midpoint
        else:
            self.lo = midpoint

    def sign(self, expression: fmpq_poly) -> int:
        """Decide an exact rational-polynomial sign at this root."""
        if expression % self.poly == 0:
            return 0
        while True:
            low, high = _interval_evaluate(expression, (self.lo, self.hi))
            if low > 0:
                return 1
            if high < 0:
                return -1
            self.refine()

    def canonical_interval(self) -> _Interval:
        """Select the first isolating centered dyadic cell."""
        if self._canonical is not None:
            return self._canonical
        variable = fmpq_poly([0, 1])
        level = 0
        while True:
            while True:
                midpoint = (self.lo + self.hi) / 2
                scaled = midpoint * (1 << level)
                index = floor(scaled + Fraction(1, 2))
                left = Fraction(2 * index - 1, 1 << (level + 1))
                right = Fraction(2 * index + 1, 1 << (level + 1))
                left_sign = self.sign(variable - _flint_polynomial((left,)))
                right_sign = self.sign(variable - _flint_polynomial((right,)))
                if left_sign == 0 or right_sign == 0:
                    break
                if left_sign > 0 and right_sign < 0:
                    if _root_count(self.sturm, left, right) == 1:
                        self._canonical = (left, right)
                        return self._canonical
                    break
                self.refine()
            level += 1


def _real_roots(P: tuple[int, ...]) -> tuple[_RealRoot, ...]:
    """Isolate and order all real roots of one irreducible polynomial."""
    polynomial = _flint_polynomial(P)
    if polynomial.degree() == 1:
        value = -Fraction(P[0], P[1])
        return (_RealRoot(P, value - 1, value + 1),)
    leading = P[-1]
    bound = 2 + max((abs(item) + leading - 1) // leading for item in P[:-1])
    chain = _sturm_chain(polynomial)
    pending = [(Fraction(-bound), Fraction(bound))]
    roots: list[_RealRoot] = []
    while pending:
        lo, hi = pending.pop()
        count = _root_count(chain, lo, hi)
        if count == 0:
            continue
        if count == 1:
            roots.append(_RealRoot(P, lo, hi))
            continue
        midpoint = (lo + hi) / 2
        if polynomial(_as_flint(midpoint)) == 0:
            raise AssertionError("irreducible nonlinear polynomial has rational root")
        pending.append((midpoint, hi))
        pending.append((lo, midpoint))
    roots.sort(key=lambda root: root.lo)
    return tuple(roots)


def _compare_distinct_roots(left: _RealRoot, right: _RealRoot) -> int:
    """Order distinct algebraic roots by refined disjoint isolators."""
    while True:
        if left.hi <= right.lo:
            return -1
        if right.hi <= left.lo:
            return 1
        if left.hi - left.lo >= right.hi - right.lo:
            left.refine()
        else:
            right.refine()


class _NumberField:
    """An exact rational number field with one selected real embedding."""

    __slots__ = ("root", "modulus", "degree")

    def __init__(self, root: _RealRoot) -> None:
        self.root = root
        self.modulus = root.poly
        self.degree = self.modulus.degree()

    def element(self, value: fmpq_poly | Fraction | int) -> _FieldElement:
        """Reduce a rational polynomial modulo this minimal polynomial."""
        polynomial = (
            value
            if isinstance(value, fmpq_poly)
            else _flint_polynomial((Fraction(value),))
        )
        return _FieldElement(self, polynomial % self.modulus)

    def theta(self) -> _FieldElement:
        """Return the selected real algebraic generator."""
        return self.element(fmpq_poly([0, 1]))


class _FieldElement:
    """An exact element of a selected real number field."""

    __slots__ = ("field", "poly")

    def __init__(self, field: _NumberField, poly: fmpq_poly) -> None:
        self.field = field
        self.poly = poly

    def _coerce(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        if isinstance(other, _FieldElement):
            if other.field is not self.field:
                raise TypeError("different number fields require an exact bridge")
            return other
        return self.field.element(other)

    def __add__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        value = self._coerce(other)
        return self.field.element(self.poly + value.poly)

    def __radd__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        return self + other

    def __neg__(self) -> _FieldElement:
        return self.field.element(-self.poly)

    def __sub__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        return self + -self._coerce(other)

    def __rsub__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        return self._coerce(other) - self

    def __mul__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        value = self._coerce(other)
        return self.field.element(self.poly * value.poly)

    def __rmul__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        return self * other

    def inverse(self) -> _FieldElement:
        """Invert a proven nonzero field element by polynomial Bezout."""
        if self.poly == 0:
            raise ZeroDivisionError("zero algebraic field element")
        common, coefficient, _ = self.poly.xgcd(self.field.modulus)
        if common.degree() != 0:
            raise AssertionError("minimal polynomial must be irreducible")
        return self.field.element(coefficient / common[0])

    def __truediv__(self, other: _FieldElement | Fraction | int) -> _FieldElement:
        return self * self._coerce(other).inverse()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, _FieldElement):
            return NotImplemented
        return self.field is other.field and self.poly == other.poly

    def is_zero(self) -> bool:
        """Decide field equality to zero without a tolerance."""
        return self.poly == 0

    def sign(self) -> int:
        """Decide the real sign under this field's selected embedding."""
        return self.field.root.sign(self.poly)

    def coefficients(self) -> _Poly:
        """Return the unique reduced ascending rational representative."""
        return _poly_coefficients(self.poly)


def _factor_minimal(polynomial: fmpq_poly) -> tuple[tuple[int, ...], ...]:
    """Factor exactly over Q and normalize each irreducible factor."""
    if polynomial == 0:
        raise ValueError("zero polynomial has no finite factorization")
    _, factors = polynomial.factor()
    return tuple(_integer_polynomial(factor) for factor, _ in factors)
