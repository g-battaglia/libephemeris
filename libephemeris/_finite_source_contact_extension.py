# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Exact quadratic field extension for singular contact candidates.

The singular stationarity line leaves a receiver quadratic over the
multiplier's rational number field. This module constructs both real roots,
including roots outside that field, and represents each in the canonical
one-root field required by the disconnected generator contract.

Provenance:
    Project-authored finite etale-algebra construction for
    ``validation/golden/specs/occultation_contact_residuals/``
    ``finite-source-contact-generator.md`` (SHA-256
    9ea40fb03e873f4149985a09010cebdd5c8427d759c830bd4998b37fa9270b7c).
    Rational polynomial factorization and matrix characteristic polynomials
    use the existing python-flint dependency. Root selection is certified
    with rational interval enclosures and exact Sturm counts.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt

from flint import fmpq_mat, fmpq_poly

from libephemeris._finite_source_contact_algebra import (
    _FieldElement,
    _Interval,
    _NumberField,
    _RealRoot,
    _as_flint,
    _as_fraction,
    _factor_minimal,
    _flint_polynomial,
    _interval_add,
    _interval_div,
    _interval_evaluate,
    _interval_mul,
    _interval_neg,
    _root_count,
    _sturm_chain,
)


def _field_evaluate(poly: fmpq_poly, value: _FieldElement) -> _FieldElement:
    """Evaluate an exact rational polynomial at one field element."""
    result = value.field.element(0)
    for coefficient in reversed(poly.coeffs()):
        result = result * value + _as_fraction(coefficient)
    return result


def _field_coordinates(value: _FieldElement) -> tuple[Fraction, ...]:
    """Return a padded rational coordinate vector in the field basis."""
    coefficients = value.coefficients()
    return coefficients + (Fraction(0),) * (value.field.degree - len(coefficients))


@dataclass(frozen=True, slots=True)
class _QuadraticElement:
    """One element `a+b*t` in a separable quadratic algebra over K."""

    algebra: _QuadraticAlgebra
    a: _FieldElement
    b: _FieldElement

    def __add__(self, other: _QuadraticElement) -> _QuadraticElement:
        if other.algebra is not self.algebra:
            raise TypeError("quadratic algebra mismatch")
        return _QuadraticElement(self.algebra, self.a + other.a, self.b + other.b)

    def __mul__(self, other: _QuadraticElement) -> _QuadraticElement:
        if other.algebra is not self.algebra:
            raise TypeError("quadratic algebra mismatch")
        return _QuadraticElement(
            self.algebra,
            self.a * other.a + self.b * other.b * self.algebra.constant,
            self.a * other.b
            + self.b * other.a
            + self.b * other.b * self.algebra.linear,
        )

    def coordinates(self) -> tuple[Fraction, ...]:
        """Return coordinates in `1,lambda,...,t,lambda*t,...` order."""
        return _field_coordinates(self.a) + _field_coordinates(self.b)


class _QuadraticAlgebra:
    """The exact K algebra satisfying `t^2=linear*t+constant`."""

    __slots__ = ("field", "linear", "constant", "dimension")

    def __init__(
        self, alpha: _FieldElement, delta: _FieldElement, epsilon: _FieldElement
    ) -> None:
        if alpha.is_zero():
            raise ValueError("receiver quadratic leading coefficient is zero")
        self.field = alpha.field
        self.linear = -2 * delta / alpha
        self.constant = -epsilon / alpha
        self.dimension = 2 * self.field.degree

    def element(
        self, a: _FieldElement | int, b: _FieldElement | int = 0
    ) -> _QuadraticElement:
        """Build one algebra element with coefficients in K."""
        first = a if isinstance(a, _FieldElement) else self.field.element(a)
        second = b if isinstance(b, _FieldElement) else self.field.element(b)
        if first.field is not self.field or second.field is not self.field:
            raise TypeError("quadratic coefficient field mismatch")
        return _QuadraticElement(self, first, second)

    def theta(self, coefficient: int) -> _QuadraticElement:
        """Construct the trial primitive element `lambda+coefficient*t`."""
        return self.element(self.field.theta(), self.field.element(coefficient))

    def basis(self) -> tuple[_QuadraticElement, ...]:
        """Return the fixed rational basis of the quadratic algebra."""
        powers = [self.field.element(1)]
        for _ in range(1, self.field.degree):
            powers.append(powers[-1] * self.field.theta())
        return tuple(self.element(power) for power in powers) + tuple(
            self.element(0, power) for power in powers
        )

    def multiplication_matrix(self, theta: _QuadraticElement) -> fmpq_mat:
        """Build exact multiplication by theta in the fixed rational basis."""
        columns = [(theta * basis).coordinates() for basis in self.basis()]
        return fmpq_mat(
            [
                [_as_flint(columns[column][row]) for column in range(self.dimension)]
                for row in range(self.dimension)
            ]
        )

    def power_matrix(self, theta: _QuadraticElement) -> fmpq_mat:
        """Use powers of a primitive element as rational basis columns."""
        power = self.element(1)
        columns = []
        for _ in range(self.dimension):
            columns.append(power.coordinates())
            power = power * theta
        return fmpq_mat(
            [
                [_as_flint(columns[column][row]) for column in range(self.dimension)]
                for row in range(self.dimension)
            ]
        )


def _sqrt_floor(value: Fraction, bits: int) -> int:
    """Floor `2**bits * sqrt(value)` by integer arithmetic."""
    if value < 0:
        raise ValueError("square-root interval must be nonnegative")
    scaled = (value.numerator << (2 * bits)) // value.denominator
    return isqrt(scaled)


def _sqrt_interval(value: _Interval, bits: int) -> _Interval:
    """Enclose square roots of an exact positive rational interval."""
    if value[0] < 0:
        raise ValueError("square-root lower bound is negative")
    denominator = 1 << bits
    lower = _sqrt_floor(value[0], bits)
    upper = _sqrt_floor(value[1], bits)
    if upper * upper * value[1].denominator < (value[1].numerator << (2 * bits)):
        upper += 1
    return Fraction(lower, denominator), Fraction(upper, denominator)


def _root_expression_interval(value: _FieldElement) -> _Interval:
    """Enclose a K element using its selected lambda isolator."""
    root = value.field.root
    return _interval_evaluate(value.poly, (root.lo, root.hi))


def _trial_theta_interval(
    alpha: _FieldElement,
    delta: _FieldElement,
    discriminant: _FieldElement,
    coefficient: int,
    branch_sign: int,
    bits: int,
) -> _Interval | None:
    """Enclose a chosen real pair's trial primitive value."""
    alpha_interval = _root_expression_interval(alpha)
    delta_interval = _root_expression_interval(delta)
    discr_interval = _root_expression_interval(discriminant)
    if alpha_interval[0] <= 0 or discr_interval[0] <= 0:
        return None
    radical = _sqrt_interval(discr_interval, bits)
    signed = radical if branch_sign > 0 else _interval_neg(radical)
    t_interval = _interval_div(
        _interval_add(_interval_neg(delta_interval), signed), alpha_interval
    )
    lambda_root = alpha.field.root
    return _interval_add(
        (lambda_root.lo, lambda_root.hi),
        _interval_mul((Fraction(coefficient),) * 2, t_interval),
    )


def _select_trial_root(
    characteristic: fmpq_poly,
    factors: tuple[tuple[int, ...], ...],
    alpha: _FieldElement,
    delta: _FieldElement,
    discriminant: _FieldElement,
    coefficient: int,
    branch_sign: int,
) -> _RealRoot:
    """Identify the exact theta root belonging to one real quadratic branch."""
    characteristic_chain = _sturm_chain(characteristic)
    bits = 8
    while True:
        interval = _trial_theta_interval(
            alpha, delta, discriminant, coefficient, branch_sign, bits
        )
        if interval is not None:
            margin = Fraction(1, 1 << bits)
            lo, hi = interval[0] - margin, interval[1] + margin
            if (
                characteristic(_as_flint(lo)) != 0
                and characteristic(_as_flint(hi)) != 0
            ):
                if _root_count(characteristic_chain, lo, hi) == 1:
                    selected = [
                        factor
                        for factor in factors
                        if _root_count(_sturm_chain(_flint_polynomial(factor)), lo, hi)
                        == 1
                    ]
                    if len(selected) != 1:
                        raise AssertionError("unique theta has no unique factor")
                    return _RealRoot(selected[0], lo, hi)
        alpha.field.root.refine()
        bits *= 2


def _global_coordinate_polynomials(
    algebra: _QuadraticAlgebra, theta: _QuadraticElement
) -> tuple[fmpq_poly, fmpq_poly]:
    """Express lambda and t as rational polynomials in a primitive theta."""
    matrix = algebra.power_matrix(theta)
    lambda_column = algebra.element(algebra.field.theta()).coordinates()
    t_column = algebra.element(0, 1).coordinates()
    right = fmpq_mat(
        [
            [_as_flint(lambda_column[row]), _as_flint(t_column[row])]
            for row in range(algebra.dimension)
        ]
    )
    solution = matrix.solve(right)
    return (
        fmpq_poly([solution[row, 0] for row in range(algebra.dimension)]),
        fmpq_poly([solution[row, 1] for row in range(algebra.dimension)]),
    )


def _selected_split_root(
    field: _NumberField,
    selected: _NumberField,
    lambda_theta: _FieldElement,
    t_theta: _FieldElement,
) -> _FieldElement:
    """Express a selected split quadratic root back in Q[lambda]."""
    degree = field.degree
    powers = [selected.element(1)]
    for _ in range(1, degree):
        powers.append(powers[-1] * lambda_theta)
    matrix = fmpq_mat(
        [
            [
                _as_flint(_field_coordinates(powers[column])[row])
                for column in range(degree)
            ]
            for row in range(degree)
        ]
    )
    right = fmpq_mat([[_as_flint(value)] for value in _field_coordinates(t_theta)])
    solution = matrix.solve(right)
    result = field.element(fmpq_poly([solution[i, 0] for i in range(degree)]))
    return result


def _quadratic_real_roots(
    alpha: _FieldElement, delta: _FieldElement, epsilon: _FieldElement
) -> tuple[tuple[_NumberField, _FieldElement, _FieldElement], ...]:
    """Construct ordered real roots and their canonical common fields.

    Each output is `(field, lambda, t)`. The field's selected generator is
    `lambda` when the root lies in Q[lambda], otherwise the first positive
    integer primitive combination `lambda+c*t`.
    """
    field = alpha.field
    if delta.field is not field or epsilon.field is not field:
        raise TypeError("quadratic coefficients must share a field")
    if alpha.sign() <= 0:
        raise ValueError("receiver quadratic requires positive leading value")
    discriminant = delta * delta - alpha * epsilon
    sign = discriminant.sign()
    if sign < 0:
        return ()
    if sign == 0:
        return ((field, field.theta(), -delta / alpha),)

    algebra = _QuadraticAlgebra(alpha, delta, epsilon)
    coefficient = 1
    while True:
        trial = algebra.theta(coefficient)
        characteristic = algebra.multiplication_matrix(trial).charpoly()
        if characteristic.gcd(characteristic.derivative()).degree() == 0:
            break
        coefficient += 1
    factors = _factor_minimal(characteristic)
    lambda_poly, t_poly = _global_coordinate_polynomials(algebra, trial)
    roots = []
    for branch_sign in (-1, 1):
        selected_root = _select_trial_root(
            characteristic,
            factors,
            alpha,
            delta,
            discriminant,
            coefficient,
            branch_sign,
        )
        selected = _NumberField(selected_root)
        lambda_theta = selected.element(lambda_poly)
        t_theta = selected.element(t_poly)
        if (
            lambda_theta + coefficient * t_theta != selected.theta()
            or not _field_evaluate(field.modulus, lambda_theta).is_zero()
        ):
            raise AssertionError("primitive coordinate reconstruction failed")
        if selected.degree == field.degree:
            t_value = _selected_split_root(field, selected, lambda_theta, t_theta)
            if (
                alpha * t_value * t_value + 2 * delta * t_value + epsilon
            ).is_zero() is False:
                raise AssertionError("split quadratic root failed its equation")
            roots.append((field, field.theta(), t_value))
        elif selected.degree == 2 * field.degree:
            roots.append((selected, lambda_theta, t_theta))
        else:
            raise AssertionError("quadratic root has unexpected field degree")
    return tuple(roots)
