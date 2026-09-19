# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Disconnected exact candidate generation for finite-source contact.

The generator exhausts the nonsingular smooth, positive-singular smooth and
core-apex branches for one exact rational target. Its decisive contact
payload must pass the independent algebraic verifier. No astronomical state,
runtime budget, public route or temporal search calls this module.

Provenance:
    Project-authored implementation of the independently reviewed private
    ``validation/golden/specs/occultation_contact_residuals/``
    ``finite-source-contact-generator.md`` (SHA-256
    9ea40fb03e873f4149985a09010cebdd5c8427d759c830bd4998b37fa9270b7c).
    Polynomial factorization and field arithmetic use the project's existing
    python-flint dependency. The verifier uses separate Fraction arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from fractions import Fraction
from functools import cmp_to_key
from itertools import combinations
from typing import Iterable, TypeAlias, TypeVar

from flint import fmpq_poly

from libephemeris._finite_source_algebraic_contact import (
    _AlgebraicContactPayload,
    _verify_algebraic_contact,
)
from libephemeris._finite_source_certificates import _Target, _ValidatedGeometry
from libephemeris._finite_source_contact_algebra import (
    _FieldElement,
    _NumberField,
    _RealRoot,
    _compare_distinct_roots,
    _factor_minimal,
    _flint_polynomial,
    _real_roots,
)
from libephemeris._finite_source_contact_extension import (
    _field_evaluate,
    _quadratic_real_roots,
)

_RationalVector: TypeAlias = tuple[Fraction, Fraction, Fraction]
_RationalMatrix: TypeAlias = tuple[_RationalVector, _RationalVector, _RationalVector]
_FieldVector: TypeAlias = tuple[_FieldElement, _FieldElement, _FieldElement]
_Poly: TypeAlias = tuple[Fraction, ...]
_ZERO = Fraction(0)
_ONE = Fraction(1)
_T = TypeVar("_T")


def _triple(values: Iterable[_T]) -> tuple[_T, _T, _T]:
    """Preserve an exact three-coordinate shape through derived expressions."""
    items = tuple(values)
    if len(items) != 3:
        raise _ContactGenerationInvariantError(
            "derived vector is not three-dimensional"
        )
    return items[0], items[1], items[2]


def _poly_const(value: Fraction | int) -> fmpq_poly:
    """Embed one exact rational in the producer's polynomial ring."""
    return _flint_polynomial((Fraction(value),))


class _ContactGenerationStatus(Enum):
    """Two private outcomes of a completed disconnected search."""

    CONTACT = auto()
    NO_CONTACT_EXHAUSTED = auto()


class _ContactBranch(Enum):
    """Exhaustive per-target contact branches in canonical order."""

    NONSINGULAR = auto()
    SINGULAR = auto()
    APEX = auto()


class _SummaryReason(Enum):
    """Why one completed branch has no point candidates."""

    NONE = auto()
    INCONSISTENT_LINEAR_SYSTEM = auto()
    NO_REAL_QUADRATIC_POINT = auto()
    PENUMBRAL_APEX_EXCLUDED = auto()
    APEX_OFF_RECEIVER = auto()


class _CandidateDisposition(Enum):
    """First exact predicate disposing of one enumerated candidate."""

    TRANSFER_SINGULAR = auto()
    FAILED_F = auto()
    FAILED_SIDE = auto()
    FAILED_GATE = auto()
    FAILED_APEX_DUAL = auto()
    ACCEPTED = auto()


class _CandidatePredicate(Enum):
    """Exact scalar recorded for a nonaccepted candidate."""

    D_ZERO = auto()
    F_NONZERO = auto()
    SIDE_NONPOSITIVE = auto()
    GATE_NEGATIVE = auto()
    TAU_NEGATIVE = auto()
    DUAL_SLACK_NEGATIVE = auto()


class _ContactGenerationInvariantError(RuntimeError):
    """The generator's derived exact identities or verifier disagree."""


@dataclass(frozen=True, slots=True)
class _ContactRoot:
    """One canonical minimal polynomial and selected real-root interval."""

    P: tuple[int, ...]
    lo: Fraction
    hi: Fraction


@dataclass(frozen=True, slots=True)
class _BranchSummary:
    """Exact counts and any branch-level absence reason."""

    branch: _ContactBranch
    completed: bool
    root_count: int
    point_count: int
    root: _ContactRoot | None
    reason: _SummaryReason
    reason_value: _Poly | None


@dataclass(frozen=True, slots=True)
class _CandidateAttempt:
    """One exact candidate or determinant-zero transfer."""

    branch: _ContactBranch
    ordinal: int
    root: _ContactRoot
    point: tuple[_Poly, _Poly, _Poly] | None
    disposition: _CandidateDisposition
    predicate: _CandidatePredicate | None
    predicate_value: _Poly | None
    payload: _AlgebraicContactPayload | None


@dataclass(frozen=True, slots=True)
class _BranchTrace:
    """A complete branch summary and every ordered candidate attempt."""

    summary: _BranchSummary
    attempts: tuple[_CandidateAttempt, ...]


@dataclass(frozen=True, slots=True)
class _ContactGenerationResult:
    """Completed private per-target contact candidate search."""

    status: _ContactGenerationStatus
    payload: _AlgebraicContactPayload | None
    trace: tuple[_BranchTrace, _BranchTrace, _BranchTrace]


@dataclass(frozen=True, slots=True)
class _TargetForms:
    """Generator-owned rational affine cone and physical gate forms."""

    a0: Fraction
    v: _RationalVector
    z0: _RationalVector
    Z: _RationalMatrix
    g0: Fraction | None
    t: _RationalVector | None
    H: Fraction


def _dot(left: _RationalVector, right: _RationalVector) -> Fraction:
    """Return an exact three-dimensional rational dot product."""
    return sum((left[i] * right[i] for i in range(3)), _ZERO)


def _cross(left: _RationalVector, right: _RationalVector) -> _RationalVector:
    """Return the right-handed rational cross product."""
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _forms(geometry: _ValidatedGeometry, target: _Target) -> _TargetForms:
    """Derive signed target forms without importing verifier helpers."""
    u = geometry.u
    w0 = -_dot(u, geometry.p)
    z0 = _triple(-item for item in _cross(u, geometry.p))
    Z: _RationalMatrix = (
        (_ZERO, -u[2], u[1]),
        (u[2], _ZERO, -u[0]),
        (-u[1], u[0], _ZERO),
    )
    s = geometry.r + geometry.m
    d_core = geometry.r - geometry.m
    if target is _Target.PENUMBRA:
        return _TargetForms(
            geometry.U * geometry.m + s * w0,
            _triple(s * item for item in u),
            z0,
            Z,
            w0 + geometry.m * s,
            u,
            geometry.U - s * s,
        )
    if target is _Target.UMBRA:
        return _TargetForms(
            geometry.U * geometry.m - d_core * w0,
            _triple(-d_core * item for item in u),
            z0,
            Z,
            w0 - geometry.m * d_core,
            u,
            geometry.U - d_core * d_core,
        )
    return _TargetForms(
        d_core * w0 - geometry.U * geometry.m,
        _triple(d_core * item for item in u),
        z0,
        Z,
        None,
        None,
        geometry.U - d_core * d_core,
    )


def _quadratic_coefficients(
    forms: _TargetForms,
) -> tuple[_RationalMatrix, _RationalVector, Fraction]:
    """Expand the signed cone residual `a^2-H*z^T*z` exactly."""
    C: _RationalMatrix = _triple(
        _triple(
            forms.v[i] * forms.v[j]
            - forms.H * sum(forms.Z[k][i] * forms.Z[k][j] for k in range(3))
            for j in range(3)
        )
        for i in range(3)
    )
    d: _RationalVector = _triple(
        forms.a0 * forms.v[i]
        - forms.H * sum(forms.Z[k][i] * forms.z0[k] for k in range(3))
        for i in range(3)
    )
    e = forms.a0 * forms.a0 - forms.H * _dot(forms.z0, forms.z0)
    return C, d, e


def _determinant(matrix: tuple[tuple[fmpq_poly, ...], ...]) -> fmpq_poly:
    """Expand a three-by-three polynomial determinant exactly."""
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _adjugate(
    matrix: tuple[tuple[fmpq_poly, ...], ...],
) -> tuple[tuple[fmpq_poly, ...], ...]:
    """Return the adjugate of a three-by-three polynomial matrix."""
    result = []
    for i in range(3):
        row = []
        for j in range(3):
            rows = [k for k in range(3) if k != j]
            columns = [k for k in range(3) if k != i]
            minor = (
                matrix[rows[0]][columns[0]] * matrix[rows[1]][columns[1]]
                - matrix[rows[0]][columns[1]] * matrix[rows[1]][columns[0]]
            )
            row.append(minor if (i + j) % 2 == 0 else -minor)
        result.append(tuple(row))
    return tuple(result)


def _stationarity_polynomials(
    geometry: _ValidatedGeometry,
    C: _RationalMatrix,
    d: _RationalVector,
    e: Fraction,
) -> tuple[fmpq_poly, tuple[fmpq_poly, ...], fmpq_poly, fmpq_poly]:
    """Construct D, N and both degree-at-most-six boundary polynomials."""
    variable = fmpq_poly([0, 1])
    M = tuple(
        tuple(
            _poly_const(geometry.A[i][j]) - variable * _poly_const(C[i][j])
            for j in range(3)
        )
        for i in range(3)
    )
    D = _determinant(M)
    adj = _adjugate(M)
    N = tuple(
        variable * sum((adj[i][j] * _poly_const(d[j]) for j in range(3)), fmpq_poly())
        for i in range(3)
    )
    Pq = (
        sum(
            (
                _poly_const(geometry.A[i][j]) * N[i] * N[j]
                for i in range(3)
                for j in range(3)
            ),
            fmpq_poly(),
        )
        - D * D
    )
    Pf = (
        sum(
            (_poly_const(C[i][j]) * N[i] * N[j] for i in range(3) for j in range(3)),
            fmpq_poly(),
        )
        + 2 * D * sum((_poly_const(d[i]) * N[i] for i in range(3)), fmpq_poly())
        + _poly_const(e) * D * D
    )
    if Pq(0) >= 0 or Pq.degree() > 6 or Pf.degree() > 6:
        raise _ContactGenerationInvariantError("stationarity polynomial bounds failed")
    return D, N, Pq, Pf


def _root_record(root: _RealRoot) -> _ContactRoot:
    """Capture the canonical selected root without a wire serialization."""
    lo, hi = root.canonical_interval()
    return _ContactRoot(root.P, lo, hi)


def _field_point(point: _FieldVector) -> tuple[_Poly, _Poly, _Poly]:
    """Return reduced coordinate polynomials for one algebraic point."""
    return _triple(item.coefficients() for item in point)


def _value(
    geometry: _ValidatedGeometry, forms: _TargetForms, point: _FieldVector
) -> tuple[
    _FieldElement, tuple[_FieldElement, ...], _FieldElement | None, _FieldElement
]:
    """Evaluate cone side, transverse vector, gate and receiver quadratic."""
    field = point[0].field
    a = field.element(forms.a0) + sum(
        (forms.v[i] * point[i] for i in range(3)), field.element(0)
    )
    z = tuple(
        field.element(forms.z0[i])
        + sum((forms.Z[i][j] * point[j] for j in range(3)), field.element(0))
        for i in range(3)
    )
    gate = None
    if forms.g0 is not None and forms.t is not None:
        gate = field.element(forms.g0) + sum(
            (forms.t[i] * point[i] for i in range(3)), field.element(0)
        )
    q = sum(
        (geometry.A[i][j] * point[i] * point[j] for i in range(3) for j in range(3)),
        field.element(0),
    )
    return a, z, gate, q


def _attempt(
    branch: _ContactBranch,
    ordinal: int,
    field: _NumberField,
    point: _FieldVector | None,
    disposition: _CandidateDisposition,
    predicate: _CandidatePredicate | None = None,
    predicate_value: _FieldElement | None = None,
    payload: _AlgebraicContactPayload | None = None,
) -> _CandidateAttempt:
    """Build one immutable exact candidate record."""
    return _CandidateAttempt(
        branch,
        ordinal,
        _root_record(field.root),
        _field_point(point) if point is not None else None,
        disposition,
        predicate,
        predicate_value.coefficients() if predicate_value is not None else None,
        payload,
    )


def _payload(
    target: _Target,
    field: _NumberField,
    tau: _FieldElement,
    beta: tuple[_FieldElement, ...],
) -> _AlgebraicContactPayload:
    """Encode derived dual values in one canonical algebraic field."""
    root = _root_record(field.root)
    return _AlgebraicContactPayload(
        root.P,
        root.lo,
        root.hi,
        tau.coefficients(),
        _triple(item.coefficients() for item in beta),
        None if target is _Target.ANTUMBRA else (_ZERO,),
    )


def _verify(
    geometry: _ValidatedGeometry,
    target: _Target,
    payload: _AlgebraicContactPayload,
) -> bool:
    """Call only the independent verifier's top-level exact checker."""
    try:
        return _verify_algebraic_contact(geometry, target, payload)
    except (TypeError, ValueError, ZeroDivisionError) as error:
        raise _ContactGenerationInvariantError(
            "generated payload is malformed"
        ) from error


def _smooth_attempt(
    geometry: _ValidatedGeometry,
    target: _Target,
    forms: _TargetForms,
    C: _RationalMatrix,
    d: _RationalVector,
    branch: _ContactBranch,
    ordinal: int,
    multiplier: _FieldElement,
    point: _FieldVector,
) -> _CandidateAttempt:
    """Apply defining checks, physical filters and independent verification."""
    field = multiplier.field
    if multiplier.sign() <= 0:
        raise _ContactGenerationInvariantError("smooth multiplier is not positive")
    a, z, gate, q = _value(geometry, forms, point)
    for i in range(3):
        left = sum((geometry.A[i][j] * point[j] for j in range(3)), field.element(0))
        right = multiplier * (
            sum((C[i][j] * point[j] for j in range(3)), field.element(d[i]))
        )
        if left != right:
            raise _ContactGenerationInvariantError(
                "smooth stationarity identity failed"
            )
    if q != field.element(1):
        raise _ContactGenerationInvariantError("smooth receiver identity failed")
    f = a * a - forms.H * sum((item * item for item in z), field.element(0))
    if f.sign() != 0:
        if branch is _ContactBranch.NONSINGULAR:
            raise _ContactGenerationInvariantError("common-root cone identity failed")
        return _attempt(
            branch,
            ordinal,
            field,
            point,
            _CandidateDisposition.FAILED_F,
            _CandidatePredicate.F_NONZERO,
            f,
        )
    if a.sign() <= 0:
        return _attempt(
            branch,
            ordinal,
            field,
            point,
            _CandidateDisposition.FAILED_SIDE,
            _CandidatePredicate.SIDE_NONPOSITIVE,
            a,
        )
    if gate is not None and gate.sign() < 0:
        return _attempt(
            branch,
            ordinal,
            field,
            point,
            _CandidateDisposition.FAILED_GATE,
            _CandidatePredicate.GATE_NEGATIVE,
            gate,
        )
    tau = 2 * multiplier * a
    beta = tuple(-2 * multiplier * forms.H * item for item in z)
    payload = _payload(target, field, tau, beta)
    if not _verify(geometry, target, payload):
        raise _ContactGenerationInvariantError("derived smooth dual failed verifier")
    return _attempt(
        branch, ordinal, field, point, _CandidateDisposition.ACCEPTED, payload=payload
    )


def _nonsingular_branch(
    geometry: _ValidatedGeometry,
    target: _Target,
    forms: _TargetForms,
    C: _RationalMatrix,
    d: _RationalVector,
    D: fmpq_poly,
    N: tuple[fmpq_poly, ...],
    Pq: fmpq_poly,
    Pf: fmpq_poly,
    singular_root: _RealRoot,
) -> _BranchTrace:
    """Enumerate every positive common root with nonsingular stationarity."""
    common = Pq.gcd(Pf)
    roots = [
        root
        for factor in _factor_minimal(common)
        for root in _real_roots(factor)
        if root.sign(fmpq_poly([0, 1])) > 0
    ]
    roots.sort(key=cmp_to_key(_compare_distinct_roots))
    attempts = []
    for ordinal, root in enumerate(roots):
        field = _NumberField(root)
        multiplier = field.theta()
        determinant = field.element(D)
        if determinant.is_zero():
            if _root_record(root) != _root_record(singular_root):
                raise _ContactGenerationInvariantError(
                    "singular transfer root mismatch"
                )
            attempts.append(
                _attempt(
                    _ContactBranch.NONSINGULAR,
                    ordinal,
                    field,
                    None,
                    _CandidateDisposition.TRANSFER_SINGULAR,
                    _CandidatePredicate.D_ZERO,
                    determinant,
                )
            )
            continue
        point: _FieldVector = _triple(field.element(item) / determinant for item in N)
        attempts.append(
            _smooth_attempt(
                geometry,
                target,
                forms,
                C,
                d,
                _ContactBranch.NONSINGULAR,
                ordinal,
                multiplier,
                point,
            )
        )
    summary = _BranchSummary(
        _ContactBranch.NONSINGULAR,
        True,
        len(roots),
        sum(item.point is not None for item in attempts),
        None,
        _SummaryReason.NONE,
        None,
    )
    return _BranchTrace(summary, tuple(attempts))


def _bridge(value: _FieldElement, multiplier: _FieldElement) -> _FieldElement:
    """Map a K element to the candidate field through lambda's polynomial."""
    if value.field is multiplier.field:
        return value
    return _field_evaluate(value.poly, multiplier)


def _singular_branch(
    geometry: _ValidatedGeometry,
    target: _Target,
    forms: _TargetForms,
    C: _RationalMatrix,
    d: _RationalVector,
    root: _RealRoot,
) -> _BranchTrace:
    """Solve the positive determinant-zero line and receiver quadratic."""
    field = _NumberField(root)
    multiplier = field.theta()
    M = tuple(
        tuple(field.element(geometry.A[i][j]) - multiplier * C[i][j] for j in range(3))
        for i in range(3)
    )
    chosen = None
    for rows in combinations(range(3), 2):
        for columns in combinations(range(3), 2):
            determinant = (
                M[rows[0]][columns[0]] * M[rows[1]][columns[1]]
                - M[rows[0]][columns[1]] * M[rows[1]][columns[0]]
            )
            if not determinant.is_zero():
                chosen = rows, columns, determinant
                break
        if chosen is not None:
            break
    if chosen is None:
        raise _ContactGenerationInvariantError(
            "positive singular matrix lacks rank two"
        )
    rows, columns, determinant = chosen
    free = next(index for index in range(3) if index not in columns)
    unused = next(index for index in range(3) if index not in rows)
    rhs = (multiplier * d[rows[0]], multiplier * d[rows[1]])
    affine = (-M[rows[0]][free], -M[rows[1]][free])

    def solve(
        vector: tuple[_FieldElement, _FieldElement],
    ) -> tuple[_FieldElement, _FieldElement]:
        return (
            (vector[0] * M[rows[1]][columns[1]] - vector[1] * M[rows[0]][columns[1]])
            / determinant,
            (M[rows[0]][columns[0]] * vector[1] - M[rows[1]][columns[0]] * vector[0])
            / determinant,
        )

    x0 = [field.element(0) for _ in range(3)]
    v0 = [field.element(0) for _ in range(3)]
    x0[columns[0]], x0[columns[1]] = solve(rhs)
    v0[columns[0]], v0[columns[1]] = solve(affine)
    v0[free] = field.element(1)
    residual_coefficient = sum(
        (M[unused][i] * v0[i] for i in range(3)), field.element(0)
    )
    if not residual_coefficient.is_zero():
        raise _ContactGenerationInvariantError(
            "singular line has nonzero residual slope"
        )
    residual = (
        sum((M[unused][i] * x0[i] for i in range(3)), field.element(0))
        - multiplier * d[unused]
    )
    if not residual.is_zero():
        return _BranchTrace(
            _BranchSummary(
                _ContactBranch.SINGULAR,
                True,
                1,
                0,
                _root_record(root),
                _SummaryReason.INCONSISTENT_LINEAR_SYSTEM,
                residual.coefficients(),
            ),
            (),
        )

    def metric(left: list[_FieldElement], right: list[_FieldElement]) -> _FieldElement:
        return sum(
            (geometry.A[i][j] * left[i] * right[j] for i in range(3) for j in range(3)),
            field.element(0),
        )

    alpha = metric(v0, v0)
    delta = metric(v0, x0)
    epsilon = metric(x0, x0) - 1
    if alpha.sign() <= 0:
        raise _ContactGenerationInvariantError(
            "singular line receiver quadratic is not convex"
        )
    discr = delta * delta - alpha * epsilon
    if discr.sign() < 0:
        return _BranchTrace(
            _BranchSummary(
                _ContactBranch.SINGULAR,
                True,
                1,
                0,
                _root_record(root),
                _SummaryReason.NO_REAL_QUADRATIC_POINT,
                discr.coefficients(),
            ),
            (),
        )
    attempts = []
    for ordinal, (candidate_field, candidate_multiplier, t_value) in enumerate(
        _quadratic_real_roots(alpha, delta, epsilon)
    ):
        point: _FieldVector = _triple(
            _bridge(x0[i], candidate_multiplier)
            + t_value * _bridge(v0[i], candidate_multiplier)
            for i in range(3)
        )
        attempts.append(
            _smooth_attempt(
                geometry,
                target,
                forms,
                C,
                d,
                _ContactBranch.SINGULAR,
                ordinal,
                candidate_multiplier,
                point,
            )
        )
    return _BranchTrace(
        _BranchSummary(
            _ContactBranch.SINGULAR,
            True,
            1,
            len(attempts),
            _root_record(root),
            _SummaryReason.NONE,
            None,
        ),
        tuple(attempts),
    )


def _apex_branch(
    geometry: _ValidatedGeometry, target: _Target, forms: _TargetForms
) -> _BranchTrace:
    """Check the rational common core apex, including opposite dual sign."""
    field = _NumberField(_RealRoot((0, 1), Fraction(-1), Fraction(1)))
    root = _root_record(field.root)
    if target is _Target.PENUMBRA:
        return _BranchTrace(
            _BranchSummary(
                _ContactBranch.APEX,
                True,
                0,
                0,
                root,
                _SummaryReason.PENUMBRAL_APEX_EXCLUDED,
                None,
            ),
            (),
        )
    d_core = geometry.r - geometry.m
    apex: _RationalVector = _triple(
        geometry.p[i] + geometry.m * geometry.u[i] / d_core for i in range(3)
    )
    q = sum(
        (geometry.A[i][j] * apex[i] * apex[j] for i in range(3) for j in range(3)),
        _ZERO,
    )
    if q != 1:
        return _BranchTrace(
            _BranchSummary(
                _ContactBranch.APEX,
                True,
                0,
                0,
                root,
                _SummaryReason.APEX_OFF_RECEIVER,
                (q - 1,),
            ),
            (),
        )
    point: _FieldVector = _triple(field.element(item) for item in apex)
    normal: _RationalVector = _triple(
        2 * sum((geometry.A[i][j] * apex[j] for j in range(3)), _ZERO) for i in range(3)
    )
    k = -d_core if target is _Target.UMBRA else d_core
    tau = _dot(geometry.u, normal) / (k * geometry.U)
    cross = _cross(geometry.u, normal)
    beta: _RationalVector = _triple(item / geometry.U for item in cross)
    derived_normal = _triple(
        k * tau * geometry.u[i] + sum(forms.Z[j][i] * beta[j] for j in range(3))
        for i in range(3)
    )
    constant = tau * forms.a0 + _dot(beta, forms.z0)
    if derived_normal != normal or constant != -2:
        raise _ContactGenerationInvariantError("apex dual construction identity failed")
    payload = _payload(
        target, field, field.element(tau), tuple(field.element(item) for item in beta)
    )
    if _verify(geometry, target, payload):
        attempt = _attempt(
            _ContactBranch.APEX,
            0,
            field,
            point,
            _CandidateDisposition.ACCEPTED,
            payload=payload,
        )
    else:
        slack = forms.H * tau * tau - _dot(beta, beta)
        if tau < 0:
            predicate = _CandidatePredicate.TAU_NEGATIVE
            value = tau
        elif slack < 0:
            predicate = _CandidatePredicate.DUAL_SLACK_NEGATIVE
            value = slack
        else:
            raise _ContactGenerationInvariantError("dual-feasible apex failed verifier")
        attempt = _attempt(
            _ContactBranch.APEX,
            0,
            field,
            point,
            _CandidateDisposition.FAILED_APEX_DUAL,
            predicate,
            field.element(value),
            payload,
        )
    return _BranchTrace(
        _BranchSummary(
            _ContactBranch.APEX,
            True,
            0,
            1,
            root,
            _SummaryReason.NONE,
            None,
        ),
        (attempt,),
    )


def _generate_algebraic_contact(
    geometry: _ValidatedGeometry, target: _Target
) -> _ContactGenerationResult:
    """Exhaust exact per-target contact candidates from rational geometry.

    Args:
        geometry: Exact private geometry; only base fields are trusted.
        target: One private signed finite-source target.

    Returns:
        A complete private contact result and exact diagnostic trace.

    Raises:
        TypeError: Entry object or target has the wrong exact type.
        ValueError: Exact geometry is malformed or outside its declared domain.
        _ContactGenerationInvariantError: A generated identity or independent
            verification conflicts with the mathematical construction.
    """
    if type(geometry) is not _ValidatedGeometry:
        raise TypeError("geometry must be an exact _ValidatedGeometry")
    if type(target) is not _Target:
        raise TypeError("target must be an exact _Target")
    fresh = _ValidatedGeometry(
        geometry.B, geometry.p, geometry.r, geometry.m, geometry.A, geometry.R
    )
    forms = _forms(fresh, target)
    C, d, e = _quadratic_coefficients(forms)
    D, N, Pq, Pf = _stationarity_polynomials(fresh, C, d, e)
    positive_singular = [
        root
        for factor in _factor_minimal(D)
        for root in _real_roots(factor)
        if root.sign(fmpq_poly([0, 1])) > 0
    ]
    if len(positive_singular) != 1:
        raise _ContactGenerationInvariantError("expected one positive determinant root")
    singular_root = positive_singular[0]
    nonsingular = _nonsingular_branch(
        fresh, target, forms, C, d, D, N, Pq, Pf, singular_root
    )
    singular = _singular_branch(fresh, target, forms, C, d, singular_root)
    apex = _apex_branch(fresh, target, forms)
    trace = (nonsingular, singular, apex)
    accepted = [
        attempt.payload
        for branch in trace
        for attempt in branch.attempts
        if attempt.disposition is _CandidateDisposition.ACCEPTED
    ]
    if len(accepted) > 1:
        raise _ContactGenerationInvariantError("one target produced multiple contacts")
    if accepted:
        payload = accepted[0]
        if payload is None:
            raise _ContactGenerationInvariantError("accepted contact lacks payload")
        return _ContactGenerationResult(
            _ContactGenerationStatus.CONTACT, payload, trace
        )
    return _ContactGenerationResult(
        _ContactGenerationStatus.NO_CONTACT_EXHAUSTED, None, trace
    )
