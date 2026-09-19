# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disconnected Decimal proposals for exact finite-source strict certificates.

Only ``_ValidatedGeometry.strict_reach`` and ``strict_miss`` decide the
physical result. Decimal iterates, averages, and support gaps are candidate
construction tools and have no proof authority.

Provenance:
    Project-authored implementation of the reviewed bounded proposer contract
    ``finite-source-strict-proposer.md`` (SHA-256
    a0300389313401ef69789a39e3ead7b490edbecef1b216ebdbbc731dbb9b516b).
    The exact affine derivation comes from the private rational geometry and
    the iteration follows the contract's projected saddle sequence. No
    astronomical source, reference implementation, or public consumer enters.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import (
    ROUND_HALF_EVEN,
    Context,
    Decimal,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    Underflow,
    localcontext,
)
from enum import Enum, auto
from fractions import Fraction
from math import isqrt
from typing import TypeAlias

from ._finite_source_certificates import (
    _DualPayload,
    _Target,
    _ValidatedGeometry,
    _quadratic,
)

_FVector: TypeAlias = tuple[Fraction, Fraction, Fraction]
_DVector: TypeAlias = tuple[Decimal, Decimal, Decimal]
_FAffine: TypeAlias = tuple[Fraction, Fraction, Fraction, Fraction]
_DAffine: TypeAlias = tuple[Decimal, Decimal, Decimal, Decimal]
_DualVector: TypeAlias = tuple[Decimal, ...]
_ZERO = Fraction(0)
_ONE = Fraction(1)
_DZERO = Decimal(0)
_DONE = Decimal(1)
_DECIMAL_EXCEPTIONS = (InvalidOperation, DivisionByZero, Overflow, Underflow)


class _StrictProposalStatus(Enum):
    """The exact verifier's accepted class or finite-policy nondecision."""

    REACH = auto()
    MISS = auto()
    UNRESOLVED = auto()


class _StrictProposalReason(Enum):
    """Why the bounded private proposal stopped."""

    ACCEPTED = auto()
    ITERATION_LIMIT = auto()
    DECIMAL_FAULT = auto()


class _StrictProposalKind(Enum):
    """Canonical attempt order within every committed checkpoint."""

    PRIMAL_CURRENT = auto()
    DUAL_CURRENT = auto()
    PRIMAL_AVERAGE = auto()
    DUAL_AVERAGE = auto()


class _StrictProposalInvariantError(RuntimeError):
    """The exact verifier or mutually exclusive certificate logic failed."""


class _DecimalFault(Exception):
    """A nonfinite or impossible Decimal support calculation occurred."""


@dataclass(frozen=True, slots=True)
class _StrictProposalPolicy:
    """Explicit iteration count and Decimal digits, without implicit defaults."""

    iterations: int
    decimal_digits: int


@dataclass(frozen=True, slots=True)
class _StrictProposalAttempt:
    """One exact candidate and its independent verifier decision."""

    iteration: int
    kind: _StrictProposalKind
    payload: _FVector | _DualPayload
    accepted: bool


@dataclass(frozen=True, slots=True)
class _StrictProposalResult:
    """Immutable finite-policy outcome and atomic checkpoint trace."""

    status: _StrictProposalStatus
    payload: _FVector | _DualPayload | None
    attempts: tuple[_StrictProposalAttempt, ...]
    iterations_completed: int
    reason: _StrictProposalReason


@dataclass(frozen=True, slots=True)
class _ExactFields:
    """Source-derived rational coefficients and norm bound."""

    a: _FAffine
    g: _FAffine | None
    z: tuple[_FAffine, _FAffine, _FAffine]
    H: Fraction
    A0: Fraction
    G0: Fraction
    L: Fraction


@dataclass(frozen=True, slots=True)
class _DecimalFields:
    """One local-context conversion of every proposal coefficient."""

    a: _DAffine
    g: _DAffine | None
    z: tuple[_DAffine, _DAffine, _DAffine]
    inverse: tuple[_DVector, _DVector, _DVector]
    sqrt_H: Decimal
    A0: Decimal
    eta: Decimal


def _affine(samples: tuple[Fraction, Fraction, Fraction, Fraction]) -> _FAffine:
    """Derive an exact intercept and three slopes from basis evaluations."""
    origin = samples[0]
    return (
        origin,
        samples[1] - origin,
        samples[2] - origin,
        samples[3] - origin,
    )


def _affine_bound(coefficients: _FAffine, radius: Fraction) -> Fraction:
    """Bound one affine scalar throughout the receiver's radius enclosure."""
    return abs(coefficients[0]) + radius * sum(
        (abs(coefficients[i]) for i in range(1, 4)), _ZERO
    )


def _exact_fields(geometry: _ValidatedGeometry, target: _Target) -> _ExactFields:
    """Extract target forms and a source-derived uniform dual gradient bound."""
    points: tuple[_FVector, _FVector, _FVector, _FVector] = (
        (_ZERO, _ZERO, _ZERO),
        (_ONE, _ZERO, _ZERO),
        (_ZERO, _ONE, _ZERO),
        (_ZERO, _ZERO, _ONE),
    )
    samples = tuple(geometry._target_values(target, point) for point in points)
    a = _affine((samples[0][0], samples[1][0], samples[2][0], samples[3][0]))
    z = tuple(
        _affine(
            (
                samples[0][3][j],
                samples[1][3][j],
                samples[2][3][j],
                samples[3][3][j],
            )
        )
        for j in range(3)
    )
    z_forms: tuple[_FAffine, _FAffine, _FAffine] = (z[0], z[1], z[2])
    gated = target is not _Target.ANTUMBRA
    g: _FAffine | None = None
    if gated:
        g0, g1, g2, g3 = (sample[1] for sample in samples)
        if g0 is None or g1 is None or g2 is None or g3 is None:
            raise _StrictProposalInvariantError("gated target lacks an affine gate")
        g = _affine((g0, g1, g2, g3))
    slope_radius = (
        geometry.r + geometry.m
        if target is _Target.PENUMBRA
        else geometry.r - geometry.m
    )
    H = geometry.U - slope_radius * slope_radius
    A0 = geometry.U * geometry.m
    G0 = geometry.m * (geometry.r + geometry.m)
    if H <= 0 or A0 <= 0 or G0 <= 0:
        raise _StrictProposalInvariantError("validated geometry lost positive scales")
    h_up = isqrt(H.numerator // H.denominator) + 1
    B_z = h_up * sum((_affine_bound(form, geometry.R) for form in z_forms), _ZERO) / A0
    if gated:
        assert g is not None
        B_a = _affine_bound(a, geometry.R) / A0
        B_g = _affine_bound(g, geometry.R) / G0
        L = max(_ONE, B_a + B_g + B_z)
    else:
        L = max(_ONE, B_z)
    return _ExactFields(a, g, z_forms, H, A0, G0, L)


def _context(digits: int) -> Context:
    """Build an isolated half-even context with exactly four enabled traps."""
    context = Context(
        prec=digits,
        rounding=ROUND_HALF_EVEN,
        Emin=-999999,
        Emax=999999,
        clamp=0,
    )
    for signal in context.traps:
        context.traps[signal] = False
    for signal in _DECIMAL_EXCEPTIONS:
        context.traps[signal] = True
    context.clear_flags()
    return context


def _finite(value: Decimal) -> Decimal:
    """Reject a nonfinite proposal intermediate without using its sign."""
    if not value.is_finite():
        raise _DecimalFault("nonfinite Decimal proposal")
    return value


def _decimal_fraction(value: Fraction) -> Decimal:
    """Convert an exact ratio by division inside the caller's local context."""
    return _finite(Decimal(value.numerator) / Decimal(value.denominator))


def _decimal_affine(coefficients: _FAffine, scale: Fraction) -> _DAffine:
    """Divide rational coefficients before Decimal conversion in order."""
    return (
        _decimal_fraction(coefficients[0] / scale),
        _decimal_fraction(coefficients[1] / scale),
        _decimal_fraction(coefficients[2] / scale),
        _decimal_fraction(coefficients[3] / scale),
    )


def _decimal_fields(
    geometry: _ValidatedGeometry, exact: _ExactFields, iterations: int
) -> _DecimalFields:
    """Convert all fixed coefficients once and compute the prescribed step."""
    sqrt_H = _finite(_decimal_fraction(exact.H).sqrt())
    a = _decimal_affine(exact.a, exact.A0)
    g = _decimal_affine(exact.g, exact.G0) if exact.g is not None else None
    z = tuple(
        tuple(_finite(_decimal_fraction(form[i] / exact.A0) * sqrt_H) for i in range(4))
        for form in exact.z
    )
    z_forms: tuple[_DAffine, _DAffine, _DAffine] = (z[0], z[1], z[2])  # type: ignore[assignment]
    inverse = tuple(
        tuple(_decimal_fraction(geometry.A_inverse[i][j]) for j in range(3))
        for i in range(3)
    )
    inverse_matrix: tuple[_DVector, _DVector, _DVector] = (
        inverse[0],
        inverse[1],
        inverse[2],
    )  # type: ignore[assignment]
    A0 = _decimal_fraction(exact.A0)
    eta = _finite(
        Decimal(3) / _finite(_decimal_fraction(exact.L) * Decimal(iterations).sqrt())
    )
    return _DecimalFields(a, g, z_forms, inverse_matrix, sqrt_H, A0, eta)


def _eval_affine(coefficients: _DAffine, point: _DVector) -> Decimal:
    """Evaluate intercept then increasing coordinate products."""
    value = coefficients[0]
    for i in range(3):
        value = _finite(value + _finite(coefficients[i + 1] * point[i]))
    return value


def _normal(fields: _DecimalFields, dual: _DualVector) -> _DVector:
    """Build the affine F slopes in the specified term order."""
    gated = fields.g is not None
    t = dual[0] if gated else _DONE
    b = dual[1:] if gated else dual
    coefficients: list[Decimal] = []
    for i in range(4):
        value = _finite(t * fields.a[i]) if gated else fields.a[i]
        for j in range(3):
            value = _finite(value + _finite(b[j] * fields.z[j][i]))
        if gated:
            assert fields.g is not None
            value = _finite(value + _finite(_finite(_DONE - t) * fields.g[i]))
        coefficients.append(value)
    return coefficients[1], coefficients[2], coefficients[3]


def _best_response(fields: _DecimalFields, dual: _DualVector) -> _DVector:
    """Maximize one Decimal affine form using A-inverse then its support."""
    normal = _normal(fields, dual)
    if all(value == _DZERO for value in normal):
        return _DZERO, _DZERO, _DZERO
    vector: list[Decimal] = []
    for i in range(3):
        component = _DZERO
        for j in range(3):
            component = _finite(component + _finite(fields.inverse[i][j] * normal[j]))
        vector.append(component)
    support_square = _DZERO
    for i in range(3):
        support_square = _finite(support_square + _finite(normal[i] * vector[i]))
    if support_square < _DZERO or support_square == _DZERO:
        raise _DecimalFault("invalid positive-definite support square")
    root = _finite(support_square.sqrt())
    return (
        _finite(vector[0] / root),
        _finite(vector[1] / root),
        _finite(vector[2] / root),
    )


def _gradient(fields: _DecimalFields, point: _DVector) -> _DualVector:
    """Evaluate the current dual gradient after the best response."""
    if fields.g is None:
        return tuple(_eval_affine(fields.z[j], point) for j in range(3))
    first = _finite(_eval_affine(fields.a, point) - _eval_affine(fields.g, point))
    z = tuple(_eval_affine(fields.z[j], point) for j in range(3))
    return (
        first,
        z[0],
        z[1],
        z[2],
    )


def _sum_components(previous: list[Decimal], current: tuple[Decimal, ...]) -> None:
    """Update a running vector sum in increasing coordinate order."""
    for i in range(len(previous)):
        previous[i] = _finite(previous[i] + current[i])


def _average(total: list[Decimal], count: int) -> tuple[Decimal, ...]:
    """Divide each accumulated component separately at a checkpoint."""
    divisor = Decimal(count)
    return tuple(_finite(total[i] / divisor) for i in range(len(total)))


def _clamp(value: Decimal, lower: Decimal, upper: Decimal) -> Decimal:
    """Clamp a finite Decimal scalar to a closed interval."""
    return lower if value < lower else upper if value > upper else value


def _projection_objective(t: Decimal, trial_t: Decimal, rho: Decimal) -> Decimal:
    """Evaluate the gated scalar projection objective in its fixed order."""
    difference = _finite(t - trial_t)
    gap = max(_DZERO, _finite(rho - t))
    return _finite(_finite(difference * difference) + _finite(gap * gap))


def _project(
    dual: _DualVector, gradient: _DualVector, fields: _DecimalFields
) -> _DualVector:
    """Project one trial dual point onto the gated or unit-ball domain."""
    trial = tuple(
        _finite(dual[i] - _finite(fields.eta * gradient[i])) for i in range(len(dual))
    )
    gated = fields.g is not None
    trial_t = trial[0] if gated else None
    trial_b = trial[1:] if gated else trial
    rho_square = _DZERO
    for i in range(3):
        rho_square = _finite(rho_square + _finite(trial_b[i] * trial_b[i]))
    rho = _finite(rho_square.sqrt())
    if not gated:
        if rho <= _DONE:
            return trial_b
        return tuple(_finite(trial_b[i] / rho) for i in range(3))
    assert trial_t is not None
    if rho == _DZERO:
        return _clamp(trial_t, _DZERO, _DONE), _DZERO, _DZERO, _DZERO
    first = _clamp(
        _finite(_finite(trial_t + rho) / Decimal(2)),
        _DZERO,
        min(rho, _DONE),
    )
    chosen = first
    first_objective = _projection_objective(first, trial_t, rho)
    if rho <= _DONE:
        second = _clamp(trial_t, rho, _DONE)
        second_objective = _projection_objective(second, trial_t, rho)
        if second_objective < first_objective or (
            second_objective == first_objective and second < first
        ):
            chosen = second
    if rho <= chosen:
        return chosen, trial_b[0], trial_b[1], trial_b[2]
    ratio = _finite(chosen / rho)
    return (
        chosen,
        _finite(ratio * trial_b[0]),
        _finite(ratio * trial_b[1]),
        _finite(ratio * trial_b[2]),
    )


def _fraction_from_decimal(value: Decimal) -> Fraction:
    """Rationalize a finite Decimal coordinate exactly, including zeros."""
    return Fraction(_finite(value))


def _primal_payload(
    point: _DVector, geometry: _ValidatedGeometry, digits: int
) -> _FVector:
    """Move a rationalized boundary point strictly inside the ellipsoid."""
    rational: _FVector = tuple(_fraction_from_decimal(point[i]) for i in range(3))  # type: ignore[assignment]
    q = _quadratic(geometry.A, rational)
    if q > _ONE:
        return rational[0] / q, rational[1] / q, rational[2] / q
    if q == _ONE:
        scale = _ONE - Fraction(1, 10**digits)
        return scale * rational[0], scale * rational[1], scale * rational[2]
    return rational


def _dual_payload(
    dual: _DualVector, exact: _ExactFields, fields: _DecimalFields
) -> _DualPayload:
    """Rationalize, clamp, and contract a proposed dual certificate."""
    if exact.g is None:
        tau = _ONE / exact.A0
        multiplier = _ZERO
        b = dual
    else:
        t = _fraction_from_decimal(dual[0])
        t = max(_ZERO, min(_ONE, t))
        tau = t / exact.A0
        multiplier = (_ONE - t) / exact.G0
        b = dual[1:]
    beta: _FVector = tuple(
        _fraction_from_decimal(_finite(_finite(fields.sqrt_H * b[i]) / fields.A0))
        for i in range(3)
    )  # type: ignore[assignment]
    Q = tau * tau * exact.H
    V = sum((beta[i] * beta[i] for i in range(3)), _ZERO)
    if V > Q:
        scale = Q / V
        beta = scale * beta[0], scale * beta[1], scale * beta[2]
    return _DualPayload(tau, beta, multiplier)


def _checked_attempt(
    geometry: _ValidatedGeometry,
    target: _Target,
    iteration: int,
    kind: _StrictProposalKind,
    payload: _FVector | _DualPayload,
) -> _StrictProposalAttempt:
    """Ask only the exact verifier to accept a well-typed candidate."""
    try:
        if kind in (
            _StrictProposalKind.PRIMAL_CURRENT,
            _StrictProposalKind.PRIMAL_AVERAGE,
        ):
            if type(payload) is not tuple:
                raise TypeError("primal candidate has the wrong type")
            accepted = geometry.strict_reach(target, payload)
        else:
            if type(payload) is not _DualPayload:
                raise TypeError("dual candidate has the wrong type")
            accepted = geometry.strict_miss(target, payload)
    except Exception as exc:
        raise _StrictProposalInvariantError(
            "exact verifier rejected a generated payload"
        ) from exc
    if type(accepted) is not bool:
        raise _StrictProposalInvariantError("exact verifier returned a non-Boolean")
    return _StrictProposalAttempt(iteration, kind, payload, accepted)


def _checkpoint(
    geometry: _ValidatedGeometry,
    target: _Target,
    iteration: int,
    digits: int,
    current_x: _DVector,
    current_y: _DualVector,
    average_x: _DVector,
    average_y: _DualVector,
    exact: _ExactFields,
    fields: _DecimalFields,
) -> tuple[_StrictProposalAttempt, ...]:
    """Build and check four provisional attempts before committing any."""
    provisional: list[_StrictProposalAttempt] = []
    provisional.append(
        _checked_attempt(
            geometry,
            target,
            iteration,
            _StrictProposalKind.PRIMAL_CURRENT,
            _primal_payload(current_x, geometry, digits),
        )
    )
    provisional.append(
        _checked_attempt(
            geometry,
            target,
            iteration,
            _StrictProposalKind.DUAL_CURRENT,
            _dual_payload(current_y, exact, fields),
        )
    )
    provisional.append(
        _checked_attempt(
            geometry,
            target,
            iteration,
            _StrictProposalKind.PRIMAL_AVERAGE,
            _primal_payload(average_x, geometry, digits),
        )
    )
    provisional.append(
        _checked_attempt(
            geometry,
            target,
            iteration,
            _StrictProposalKind.DUAL_AVERAGE,
            _dual_payload(average_y, exact, fields),
        )
    )
    return tuple(provisional)


def _unresolved(
    attempts: list[_StrictProposalAttempt],
    completed: int,
    reason: _StrictProposalReason,
) -> _StrictProposalResult:
    """Freeze a nondecision without manufacturing a physical class."""
    return _StrictProposalResult(
        _StrictProposalStatus.UNRESOLVED, None, tuple(attempts), completed, reason
    )


def _run(
    geometry: _ValidatedGeometry, target: _Target, policy: _StrictProposalPolicy
) -> _StrictProposalResult:
    """Execute the fixed finite-precision sequence with atomic checkpoints."""
    exact = _exact_fields(geometry, target)
    attempts: list[_StrictProposalAttempt] = []
    completed = 0
    context = _context(policy.decimal_digits)
    with localcontext(context) as active:
        active.clear_flags()
        try:
            fields = _decimal_fields(geometry, exact, policy.iterations)
            dual: _DualVector = (
                (Decimal(1) / Decimal(2), _DZERO, _DZERO, _DZERO)
                if exact.g is not None
                else (_DZERO, _DZERO, _DZERO)
            )
            sum_x = [_DZERO, _DZERO, _DZERO]
            sum_y = [_DZERO for _ in range(len(dual))]
        except (*_DECIMAL_EXCEPTIONS, _DecimalFault):
            return _unresolved(attempts, completed, _StrictProposalReason.DECIMAL_FAULT)
        for iteration in range(1, policy.iterations + 1):
            try:
                point = _best_response(fields, dual)
                _sum_components(sum_x, point)
                _sum_components(sum_y, dual)
                gradient = _gradient(fields, point)
                checkpoint = (
                    iteration & (iteration - 1) == 0 or iteration == policy.iterations
                )
                provisional: tuple[_StrictProposalAttempt, ...] | None = None
                if checkpoint:
                    average_x = _average(sum_x, iteration)
                    average_y = _average(sum_y, iteration)
                    primal_average: _DVector = (
                        average_x[0],
                        average_x[1],
                        average_x[2],
                    )
                    provisional = _checkpoint(
                        geometry,
                        target,
                        iteration,
                        policy.decimal_digits,
                        point,
                        dual,
                        primal_average,
                        average_y,
                        exact,
                        fields,
                    )
            except (*_DECIMAL_EXCEPTIONS, _DecimalFault):
                return _unresolved(
                    attempts, completed, _StrictProposalReason.DECIMAL_FAULT
                )
            if provisional is not None:
                reaches = [
                    item
                    for item in provisional
                    if item.accepted
                    and item.kind
                    in (
                        _StrictProposalKind.PRIMAL_CURRENT,
                        _StrictProposalKind.PRIMAL_AVERAGE,
                    )
                ]
                misses = [
                    item
                    for item in provisional
                    if item.accepted
                    and item.kind
                    in (
                        _StrictProposalKind.DUAL_CURRENT,
                        _StrictProposalKind.DUAL_AVERAGE,
                    )
                ]
                if reaches and misses:
                    raise _StrictProposalInvariantError(
                        "both strict classes verified at one checkpoint"
                    )
                attempts.extend(provisional)
                completed = iteration
                if reaches or misses:
                    accepted = next(item for item in provisional if item.accepted)
                    status = (
                        _StrictProposalStatus.REACH
                        if reaches
                        else _StrictProposalStatus.MISS
                    )
                    return _StrictProposalResult(
                        status,
                        accepted.payload,
                        tuple(attempts),
                        completed,
                        _StrictProposalReason.ACCEPTED,
                    )
            else:
                completed = iteration
            if iteration == policy.iterations:
                return _unresolved(
                    attempts, completed, _StrictProposalReason.ITERATION_LIMIT
                )
            try:
                dual = _project(dual, gradient, fields)
            except (*_DECIMAL_EXCEPTIONS, _DecimalFault):
                return _unresolved(
                    attempts, completed, _StrictProposalReason.DECIMAL_FAULT
                )
    raise _StrictProposalInvariantError("fixed proposal loop ended without an outcome")


def _propose_strict_certificate(
    geometry: _ValidatedGeometry, target: _Target, policy: _StrictProposalPolicy
) -> _StrictProposalResult:
    """Propose one per-target certificate through a fresh exact geometry."""
    if type(target) is not _Target:
        raise TypeError("target must be an exact private _Target")
    if type(policy) is not _StrictProposalPolicy:
        raise TypeError("policy must be an exact _StrictProposalPolicy")
    if type(policy.iterations) is not int:
        raise TypeError("iterations must be an exact int")
    if type(policy.decimal_digits) is not int:
        raise TypeError("decimal_digits must be an exact int")
    if policy.iterations < 1:
        raise ValueError("iterations must be at least one")
    if policy.decimal_digits < 2:
        raise ValueError("decimal_digits must be at least two")
    if type(geometry) is not _ValidatedGeometry:
        raise TypeError("geometry must be an exact _ValidatedGeometry")
    fresh = _ValidatedGeometry(
        geometry.B, geometry.p, geometry.r, geometry.m, geometry.A, geometry.R
    )
    return _run(fresh, target, policy)
