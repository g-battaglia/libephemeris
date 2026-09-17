# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Private certified frames for cone and ellipsoid contact geometry.

A circular shadow cone needs more than a scalar axis offset when the receiving
body is an arbitrarily oriented ellipsoid. This module carries the ellipsoid
metric and the complete directed axis line as immutable source quantities. The
unit axis direction and its closest anchor are derived in Arb from one raw span
and one point on the line, rather than accepted as independently rounded floats
that merely look normalized or orthogonal.

The frame is private: it is an owned producer/consumer value, not public API and
not process-global metadata. This module supplies only input validation and the
axis-line derivation. It does not claim to solve the cone reach, KKT active sets,
nappe boundaries, or time-dependent contact problem.

Provenance:
    Project-authored representation of the coordinate-free line and ellipsoid
    definitions in standard quadric geometry. The ellipsoid ``x.T A x <= 1``
    and closest-point line decomposition are elementary Euclidean identities.
    Positive-definiteness and binary64 ingestion use the project's reviewed Arb
    interval primitives. No astronomical coefficient, compatibility output, or
    fitted threshold enters this module.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Final

from flint import arb, arb_mat

from .intervals import (
    Ball,
    BallMatrix,
    IntervalCertificationError,
    ball_from_float,
    interval_cholesky,
    krawczyk_image,
    strictly_contains_vector,
)

__all__: list[str] = []

_VECTOR_DIMENSION: Final = 3

Vector3 = tuple[float, float, float]
Matrix3 = tuple[Vector3, Vector3, Vector3]
BallVector3 = tuple[Ball, Ball, Ball]


@dataclass(frozen=True, slots=True)
class _ConeSection:
    """One selected cone section on the fundamental plane."""

    radius_km: float
    cosine: float
    branch_sign: int

    def validate(self) -> None:
        """Validate a non-negative section and acute cone half-angle."""
        if type(self.radius_km) is not float or not math.isfinite(self.radius_km):
            raise ValueError("cone radius must be a finite native float")
        if self.radius_km < 0.0:
            raise ValueError("cone radius must be non-negative")
        if type(self.cosine) is not float or not math.isfinite(self.cosine):
            raise ValueError("cone cosine must be a finite native float")
        if not 0.0 < self.cosine <= 1.0:
            raise ValueError("cone cosine must lie in (0, 1]")
        if type(self.branch_sign) is not int or self.branch_sign not in {-1, 1}:
            raise ValueError("cone branch sign must be exactly -1 or +1")


class _CoreShadowClass(str, Enum):
    """Classification carried separately from non-negative cone geometry."""

    UMBRA = "umbra"
    ANTUMBRA = "antumbra"
    APEX = "apex"


@dataclass(frozen=True, slots=True)
class _CoreConeSection:
    """Physical core cone paired with its signed-diameter classification."""

    cone: _ConeSection
    signed_diameter_km: float
    shadow_class: _CoreShadowClass


@dataclass(frozen=True, slots=True)
class _RegularContactEvaluation:
    """Regular cone/ellipsoid residuals and stationarity vector."""

    ellipsoid: Ball
    cone: Ball
    nappe_radius: Ball
    radial_distance: Ball
    stationarity: BallVector3


@dataclass(frozen=True, slots=True)
class _NappeBoundaryEvaluation:
    """KKT residuals on the active ``r=0`` boundary with ``rho>0``."""

    ellipsoid: Ball
    nappe_radius: Ball
    radial_distance: Ball
    stationarity: BallVector3


@dataclass(frozen=True, slots=True)
class _RadialSubgradientEvaluation:
    """KKT residuals for one supplied ``rho=0`` subgradient witness."""

    ellipsoid: Ball
    radial_squared: Ball
    nappe_radius: Ball
    witness_axis_dot: Ball
    witness_norm_squared: Ball
    stationarity: BallVector3


@dataclass(frozen=True, slots=True)
class _ApexEvaluation:
    """Exact geometric apex candidate for a nonzero-angle cone."""

    axial_parameter: Ball
    point_km: BallVector3
    ellipsoid: Ball
    cone: Ball


@dataclass(frozen=True, slots=True)
class _RegularRootCertificate:
    """Strict Krawczyk inclusion for one regular KKT root."""

    frame: _EllipsoidContactFrame
    cone: _ConeSection
    root_box: tuple[Ball, Ball, Ball, Ball]
    image: tuple[Ball, Ball, Ball, Ball]
    cone_residual: Ball


@dataclass(frozen=True, slots=True)
class _ResidualReduction:
    """Ordered minimum and normalized sign within a supplied candidate set."""

    minimum: Ball
    normalized: Ball
    sign: int
    candidate_count: int


class _GlobalReachStatus(str, Enum):
    """Private result classes for standalone global certification."""

    REACH = "reach"
    MISS = "miss"
    CONTACT = "contact"
    UNRESOLVED = "unresolved"
    INVALID = "invalid"


@dataclass(frozen=True, slots=True)
class _ProjectedDualWitness:
    """Exact rational projected payload for one transverse dual vector."""

    raw_span: tuple[Fraction, Fraction, Fraction]
    seed: tuple[Fraction, Fraction, Fraction]
    scale: Fraction
    multiplier: Fraction


@dataclass(frozen=True, slots=True)
class _CrossDualWitness:
    """Exact rational cross-product payload for one transverse dual vector."""

    raw_span: tuple[Fraction, Fraction, Fraction]
    seed: tuple[Fraction, Fraction, Fraction]
    scale: Fraction
    multiplier: Fraction


DualWitness = _ProjectedDualWitness | _CrossDualWitness


class _ContactEqualityVariant(str, Enum):
    """Reviewed exact-identity routes for one contact equality."""

    ZERO_ANGLE_LINE = "zero-angle-line"


class _ZeroAngleLineRelation(str, Enum):
    """Exact source-level relation between one axis and the ellipsoid."""

    CROSSING = "crossing"
    TANGENT = "tangent"


@dataclass(frozen=True, slots=True)
class _ContactEqualityProof:
    """Typed exact-equality payload independent of ordinary interval brackets."""

    variant: _ContactEqualityVariant
    dual_witness: DualWitness
    line_relation: _ZeroAngleLineRelation


@dataclass(frozen=True, slots=True)
class _GlobalReachProof:
    """Standalone primal/dual bracket for one frame and cone."""

    status: _GlobalReachStatus
    lower_bound: Ball | None
    upper_bound: Ball | None
    primal_point: BallVector3 | None
    dual_witness: DualWitness | None
    equality_proof: _ContactEqualityProof | None
    reason: str


@dataclass(frozen=True, slots=True)
class _CertifiedAxisLine:
    """Arb enclosure of a unit direction and closest axis anchor."""

    direction: BallVector3
    anchor_km: BallVector3


@dataclass(frozen=True, slots=True)
class _EllipsoidContactFrame:
    """Immutable source frame for a private ellipsoidal contact calculation.

    ``axis_point_km`` is any point on the axis. ``axis_span`` is a nonzero
    vector directed downstream along the shadow. They are retained as exact
    binary64 inputs; :meth:`certified_axis_line` defines the unit direction and
    closest anchor from them in one expression graph. Branch signs select the
    physical nappe independently for the penumbral and core cones.
    """

    metric_km_minus_2: Matrix3
    axis_point_km: Vector3
    axis_span: Vector3
    penumbra_branch_sign: int
    core_branch_sign: int

    def validate(self) -> None:
        """Validate finite native inputs and certified frame invariants.

        Raises:
            ValueError: If a scalar has the wrong type, is non-finite, a branch
                sign is not exactly ``-1`` or ``+1``, or the metric is not
                exactly symmetric.
            IntervalCertificationError: If the metric is not certified positive
                definite or the axis span is not certified nonzero.
        """
        for name, vector in (
            ("axis point", self.axis_point_km),
            ("axis span", self.axis_span),
        ):
            _validate_vector(vector, name)
        _validate_metric(self.metric_km_minus_2)
        for name, sign in (
            ("penumbra branch", self.penumbra_branch_sign),
            ("core branch", self.core_branch_sign),
        ):
            if type(sign) is not int or sign not in {-1, 1}:
                raise ValueError(f"{name} sign must be exactly -1 or +1")
        span = _ball_vector(self.axis_span)
        squared_norm = sum((component * component for component in span), arb(0))
        if not squared_norm > 0:
            raise IntervalCertificationError(
                "axis span has no certified positive squared norm"
            )

    def metric_ball_matrix(self) -> BallMatrix:
        """Return the exact-input ellipsoid metric as an Arb matrix."""
        self.validate()
        return arb_mat(
            [
                [ball_from_float(value) for value in row]
                for row in self.metric_km_minus_2
            ]
        )

    def certified_axis_line(self) -> _CertifiedAxisLine:
        """Derive the unit direction and closest anchor in Arb.

        The definitions are

        ``e = span / ||span||`` and ``q = point - (point dot e) e``.

        They make unit length and ``q dot e = 0`` construction identities rather
        than tolerance checks on separately rounded public values.
        """
        self.validate()
        span = _ball_vector(self.axis_span)
        point = _ball_vector(self.axis_point_km)
        length = sum((component * component for component in span), arb(0)).sqrt()
        if not length > 0:
            raise IntervalCertificationError("axis span length is not positive")
        direction = tuple(component / length for component in span)
        along = sum(
            (point[index] * direction[index] for index in range(_VECTOR_DIMENSION)),
            arb(0),
        )
        anchor = tuple(
            point[index] - along * direction[index]
            for index in range(_VECTOR_DIMENSION)
        )
        if any(not value.is_finite() for value in (*direction, *anchor)):
            raise IntervalCertificationError("axis-line derivation is not finite")
        return _CertifiedAxisLine(direction, anchor)  # type: ignore[arg-type]


def _core_cone_section(
    signed_diameter_km: float,
    cosine: float,
) -> _CoreConeSection:
    """Keep signed core classification separate from physical cone radius."""
    if type(signed_diameter_km) is not float or not math.isfinite(signed_diameter_km):
        raise ValueError("signed core diameter must be a finite native float")
    if signed_diameter_km < 0.0:
        shadow_class = _CoreShadowClass.UMBRA
        branch_sign = -1
    elif signed_diameter_km > 0.0:
        shadow_class = _CoreShadowClass.ANTUMBRA
        branch_sign = 1
    else:
        shadow_class = _CoreShadowClass.APEX
        branch_sign = 1
    cone = _ConeSection(abs(signed_diameter_km) / 2.0, cosine, branch_sign)
    cone.validate()
    return _CoreConeSection(cone, signed_diameter_km, shadow_class)


def _dot(left: tuple[Ball, ...], right: tuple[Ball, ...]) -> Ball:
    """Return one outward-rounded Euclidean dot product."""
    return sum((left[index] * right[index] for index in range(len(left))), arb(0))


def _metric_vector(metric: BallMatrix, vector: BallVector3) -> BallVector3:
    """Multiply one three-dimensional Arb matrix/vector pair."""
    column = arb_mat([[value] for value in vector])
    product = metric * column
    return tuple(product[index, 0] for index in range(_VECTOR_DIMENSION))  # type: ignore[return-value]


def _evaluate_smooth_contact(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    ellipsoid_multiplier: Ball,
    *,
    require_positive_nappe: bool,
) -> _RegularContactEvaluation:
    """Evaluate smooth cone geometry, optionally requiring inactive nappe.

    The radial domain ``rho > 0`` is mandatory because the cone gradient is
    otherwise nonsmooth. The nappe radius is returned as an equation; regular
    callers require it positive, while the active-boundary evaluator permits an
    enclosure containing zero. This helper evaluates no global classification.

    Args:
        frame: Validated private ellipsoid and axis source frame.
        cone: One validated physical cone section.
        point_km: Three-dimensional interval point.
        ellipsoid_multiplier: Non-negative KKT multiplier enclosure.
        require_positive_nappe: Require ``r>0`` for the inactive constraint.

    Returns:
        Outward enclosures of the ellipsoid equation, cone equation, nappe
        radius, radial distance, and three stationarity components.

    Raises:
        ValueError: If a point/multiplier/cone input is malformed or non-finite.
        IntervalCertificationError: If ``rho>0`` or a requested ``r>0`` is not
            certified.
    """
    frame.validate()
    cone.validate()
    if (
        type(point_km) is not tuple
        or len(point_km) != _VECTOR_DIMENSION
        or any(type(value) is not arb or not value.is_finite() for value in point_km)
    ):
        raise ValueError("contact point must contain three finite Arb intervals")
    if (
        type(ellipsoid_multiplier) is not arb
        or not ellipsoid_multiplier.is_finite()
        or ellipsoid_multiplier.lower() < 0
    ):
        raise ValueError("ellipsoid multiplier must be a finite non-negative interval")

    metric = frame.metric_ball_matrix()
    axis = frame.certified_axis_line()
    relative = tuple(
        point_km[index] - axis.anchor_km[index] for index in range(_VECTOR_DIMENSION)
    )
    axial = _dot(point_km, axis.direction)
    perpendicular = tuple(
        relative[index] - _dot(relative, axis.direction) * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    radial_squared = _dot(perpendicular, perpendicular)
    if not radial_squared > 0:
        raise IntervalCertificationError(
            "regular contact requires radial distance separated above zero"
        )
    radial = radial_squared.sqrt()
    cosine = ball_from_float(cone.cosine)
    sine = (1 - cosine * cosine).sqrt()
    radius = ball_from_float(cone.radius_km)
    branch = arb(cone.branch_sign)
    tangent = sine / cosine
    nappe_radius = radius + branch * tangent * axial
    if require_positive_nappe and not nappe_radius > 0:
        raise IntervalCertificationError(
            "regular contact requires nappe radius separated above zero"
        )
    metric_point = _metric_vector(metric, point_km)
    ellipsoid = _dot(point_km, metric_point) - 1
    cone_value = cosine * radial - branch * sine * axial - cosine * radius
    normal = tuple(
        cosine * perpendicular[index] / radial - branch * sine * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    stationarity = tuple(
        normal[index] + 2 * ellipsoid_multiplier * metric_point[index]
        for index in range(_VECTOR_DIMENSION)
    )
    return _RegularContactEvaluation(
        ellipsoid,
        cone_value,
        nappe_radius,
        radial,
        stationarity,  # type: ignore[arg-type]
    )


def _evaluate_regular_contact(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    ellipsoid_multiplier: Ball,
) -> _RegularContactEvaluation:
    """Evaluate the smooth KKT terms with inactive nappe constraint."""
    return _evaluate_smooth_contact(
        frame,
        cone,
        point_km,
        ellipsoid_multiplier,
        require_positive_nappe=True,
    )


def _evaluate_nappe_boundary(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    ellipsoid_multiplier: Ball,
    nappe_multiplier: Ball,
) -> _NappeBoundaryEvaluation:
    """Evaluate the active-nappe KKT system for ``r=0`` and ``rho>0``.

    The stationarity equation is

    ``grad(g) + 2*lambda*A*x - mu*k*tan(f)*e = 0``

    with both multipliers non-negative. This evaluator includes the ellipsoid
    equation and the active ``r=0`` equation; it does not isolate a root or
    prove that this active set is globally minimizing.

    Args:
        frame: Validated private ellipsoid and axis frame.
        cone: Validated cone section.
        point_km: Three finite Arb coordinate intervals.
        ellipsoid_multiplier: Certified non-negative ``lambda`` enclosure.
        nappe_multiplier: Certified non-negative ``mu`` enclosure.

    Returns:
        Outward enclosures of ``E``, ``r``, ``rho``, and stationarity.

    Raises:
        ValueError: If inputs are malformed or multipliers are not certified
            non-negative.
        IntervalCertificationError: If ``rho>0`` is not certified.
    """
    regular = _evaluate_smooth_contact(
        frame,
        cone,
        point_km,
        ellipsoid_multiplier,
        require_positive_nappe=False,
    )
    if (
        type(nappe_multiplier) is not arb
        or not nappe_multiplier.is_finite()
        or nappe_multiplier.lower() < 0
    ):
        raise ValueError("nappe multiplier must be a finite non-negative interval")
    axis = frame.certified_axis_line()
    cosine = ball_from_float(cone.cosine)
    tangent = (1 - cosine * cosine).sqrt() / cosine
    branch = arb(cone.branch_sign)
    stationarity = tuple(
        regular.stationarity[index]
        - nappe_multiplier * branch * tangent * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    return _NappeBoundaryEvaluation(
        regular.ellipsoid,
        regular.nappe_radius,
        regular.radial_distance,
        stationarity,  # type: ignore[arg-type]
    )


def _nappe_boundary_jacobian(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    ellipsoid_multiplier: Ball,
    nappe_multiplier: Ball,
) -> BallMatrix:
    """Return the Jacobian of ``[E, r, active stationarity]``.

    The unknown vector is ``(x,y,z,lambda,mu)``. The active constraint row is
    ``grad(r)=k*tan(f)*e``; the spatial stationarity block shares the regular cone
    Hessian plus ``2*lambda*A``; multiplier columns are ``2*A*x`` and
    ``-k*tan(f)*e``.
    """
    evaluation = _evaluate_nappe_boundary(
        frame, cone, point_km, ellipsoid_multiplier, nappe_multiplier
    )
    metric = frame.metric_ball_matrix()
    axis = frame.certified_axis_line()
    relative = tuple(
        point_km[index] - axis.anchor_km[index] for index in range(_VECTOR_DIMENSION)
    )
    projected = _dot(relative, axis.direction)
    perpendicular = tuple(
        relative[index] - projected * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    radial = evaluation.radial_distance
    radial_cubed = radial * radial * radial
    cosine = ball_from_float(cone.cosine)
    tangent = (1 - cosine * cosine).sqrt() / cosine
    branch = arb(cone.branch_sign)
    metric_point = _metric_vector(metric, point_km)
    jacobian = arb_mat(5, 5)
    for column in range(_VECTOR_DIMENSION):
        jacobian[0, column] = 2 * metric_point[column]
        jacobian[1, column] = branch * tangent * axis.direction[column]
    for row in range(_VECTOR_DIMENSION):
        for column in range(_VECTOR_DIMENSION):
            projector = arb(1 if row == column else 0) - (
                axis.direction[row] * axis.direction[column]
            )
            hessian = cosine * (
                projector / radial
                - perpendicular[row] * perpendicular[column] / radial_cubed
            )
            jacobian[row + 2, column] = (
                hessian + 2 * ellipsoid_multiplier * metric[row, column]
            )
        jacobian[row + 2, 3] = 2 * metric_point[row]
        jacobian[row + 2, 4] = -branch * tangent * axis.direction[row]
    return jacobian


def _evaluate_radial_subgradient(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    witness: BallVector3,
    ellipsoid_multiplier: Ball,
    nappe_multiplier: Ball,
) -> _RadialSubgradientEvaluation:
    """Evaluate KKT residuals at ``rho=0`` using a supplied subgradient.

    The supplied witness represents ``v`` in equation (7). Its admissible set is
    ``v dot e = 0`` and ``||v|| <= 1``. This function returns those constraint
    residuals together with ``rho^2``, ``E``, ``r``, and stationarity; it does
    not claim to prove the existential witness or isolate a nonsmooth root.

    Args:
        frame: Validated private ellipsoid and axis frame.
        cone: Validated selected cone section.
        point_km: Three finite Arb point intervals.
        witness: Three finite Arb subgradient intervals.
        ellipsoid_multiplier: Certified non-negative ``lambda`` enclosure.
        nappe_multiplier: Certified non-negative ``mu`` enclosure.

    Returns:
        Outward enclosures of all nonsmooth KKT residuals and witness bounds.

    Raises:
        ValueError: If inputs or multipliers are malformed or not certified
            non-negative.
    """
    frame.validate()
    cone.validate()
    for values, name in ((point_km, "point"), (witness, "subgradient witness")):
        if (
            type(values) is not tuple
            or len(values) != _VECTOR_DIMENSION
            or any(type(value) is not arb or not value.is_finite() for value in values)
        ):
            raise ValueError(f"radial {name} must contain three finite Arb intervals")
    for value, name in (
        (ellipsoid_multiplier, "ellipsoid multiplier"),
        (nappe_multiplier, "nappe multiplier"),
    ):
        if type(value) is not arb or not value.is_finite() or value.lower() < 0:
            raise ValueError(f"radial {name} must be a finite non-negative interval")

    metric = frame.metric_ball_matrix()
    axis = frame.certified_axis_line()
    relative = tuple(
        point_km[index] - axis.anchor_km[index] for index in range(_VECTOR_DIMENSION)
    )
    projected = _dot(relative, axis.direction)
    perpendicular = tuple(
        relative[index] - projected * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    radial_squared = _dot(perpendicular, perpendicular)
    axial = _dot(point_km, axis.direction)
    cosine = ball_from_float(cone.cosine)
    sine = (1 - cosine * cosine).sqrt()
    tangent = sine / cosine
    branch = arb(cone.branch_sign)
    radius = ball_from_float(cone.radius_km)
    nappe_radius = radius + branch * tangent * axial
    metric_point = _metric_vector(metric, point_km)
    ellipsoid = _dot(point_km, metric_point) - 1
    witness_axis_dot = _dot(witness, axis.direction)
    witness_norm_squared = _dot(witness, witness)
    stationarity = tuple(
        cosine * witness[index]
        - branch * sine * axis.direction[index]
        + 2 * ellipsoid_multiplier * metric_point[index]
        - nappe_multiplier * branch * tangent * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    return _RadialSubgradientEvaluation(
        ellipsoid,
        radial_squared,
        nappe_radius,
        witness_axis_dot,
        witness_norm_squared,
        stationarity,  # type: ignore[arg-type]
    )


def _evaluate_cone_apex(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
) -> _ApexEvaluation:
    """Evaluate the unique joint ``rho=0, r=0`` nonzero-angle candidate.

    Since the certified axis anchor is perpendicular to ``e``, every radial-axis
    point is ``q+t*e`` and has axial coordinate ``t``. For ``sin(f)>0``, the
    active equation gives ``t=-l/(k*tan(f))``. Since the point is constructed
    on the axis and the active radius is zero, ``rho=0`` and ``g=0`` are exact
    construction identities. No stationarity or global-minimum claim is made.

    Raises:
        ValueError: If the cone has zero angle and therefore no finite apex.
        IntervalCertificationError: If the sine is not finite and certified
            strictly positive at the current precision.
    """
    frame.validate()
    cone.validate()
    axis = frame.certified_axis_line()
    cosine = ball_from_float(cone.cosine)
    sine_squared = 1 - cosine * cosine
    if sine_squared.is_exact() and sine_squared.is_zero():
        raise ValueError("a zero-angle cone has no finite apex")
    if not sine_squared.is_finite() or not sine_squared > 0:
        raise IntervalCertificationError(
            "cone sine is not finite and separated above zero"
        )
    sine = sine_squared.sqrt()
    tangent = sine / cosine
    branch = arb(cone.branch_sign)
    radius = ball_from_float(cone.radius_km)
    axial = -radius / (branch * tangent)
    point = tuple(
        axis.anchor_km[index] + axial * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    metric = frame.metric_ball_matrix()
    metric_point = _metric_vector(metric, point)  # type: ignore[arg-type]
    ellipsoid = _dot(point, metric_point) - 1
    cone_residual = arb(0)
    return _ApexEvaluation(
        axial,
        point,  # type: ignore[arg-type]
        ellipsoid,
        cone_residual,
    )


def _regular_contact_jacobian(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    ellipsoid_multiplier: Ball,
) -> BallMatrix:
    """Return the analytic Jacobian of ``[E, stationarity]``.

    The unknown vector is ``(x, y, z, lambda)``. On the regular domain,

    ``H_g = c * (P/rho - w*w.T/rho^3)``

    and the stationarity block is ``H_g + 2*lambda*A`` with final column
    ``2*A*x``. The first row is the ellipsoid gradient ``2*A*x`` and zero in
    the multiplier column. The cone residual is evaluated separately after a
    KKT root has been isolated; it is not a fifth unknown equation.

    Args:
        frame: Validated private ellipsoid and axis source frame.
        cone: One validated physical cone section.
        point_km: Three-dimensional interval point in the regular domain.
        ellipsoid_multiplier: Certified non-negative multiplier enclosure.

    Returns:
        Outward-rounded ``4 x 4`` interval Jacobian.

    Raises:
        ValueError: If an input violates the regular evaluator contract.
        IntervalCertificationError: If the smooth-domain conditions are not
            certified.
    """
    evaluation = _evaluate_regular_contact(frame, cone, point_km, ellipsoid_multiplier)
    metric = frame.metric_ball_matrix()
    axis = frame.certified_axis_line()
    relative = tuple(
        point_km[index] - axis.anchor_km[index] for index in range(_VECTOR_DIMENSION)
    )
    projected = _dot(relative, axis.direction)
    perpendicular = tuple(
        relative[index] - projected * axis.direction[index]
        for index in range(_VECTOR_DIMENSION)
    )
    radial = evaluation.radial_distance
    radial_cubed = radial * radial * radial
    cosine = ball_from_float(cone.cosine)
    metric_point = _metric_vector(metric, point_km)

    jacobian = arb_mat(_VECTOR_DIMENSION + 1, _VECTOR_DIMENSION + 1)
    for column in range(_VECTOR_DIMENSION):
        jacobian[0, column] = 2 * metric_point[column]
    for row in range(_VECTOR_DIMENSION):
        for column in range(_VECTOR_DIMENSION):
            projector = arb(1 if row == column else 0) - (
                axis.direction[row] * axis.direction[column]
            )
            cone_hessian = cosine * (
                projector / radial
                - perpendicular[row] * perpendicular[column] / radial_cubed
            )
            jacobian[row + 1, column] = (
                cone_hessian + 2 * ellipsoid_multiplier * metric[row, column]
            )
        jacobian[row + 1, _VECTOR_DIMENSION] = 2 * metric_point[row]
    return jacobian


def _certify_regular_root_box(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    root_box: tuple[Ball, Ball, Ball, Ball],
) -> _RegularRootCertificate:
    """Certify one regular KKT root by strict Krawczyk inclusion.

    The first three variables are a point on the ellipsoid and the fourth is the
    non-negative ellipsoid multiplier. The box must already exclude ``rho=0``
    and ``r=0``. This function proves one root of ``[E, stationarity]`` inside
    the supplied box. It returns the cone residual for a separate exact-contact
    decision; it does not claim that the root lies on ``g=0`` or is the global
    minimizer.

    Args:
        frame: Validated private contact frame.
        cone: Validated selected cone section.
        root_box: Four finite Arb intervals ``(x,y,z,lambda)``.

    Returns:
        The root box, strict Krawczyk image, and cone residual at the box.

    Raises:
        ValueError: If the box shape or entries are invalid.
        IntervalCertificationError: If the midpoint Jacobian is singular, the
            Krawczyk image is not strictly inside the box, or the smooth domain
            is not certified.
    """
    if (
        type(root_box) is not tuple
        or len(root_box) != _VECTOR_DIMENSION + 1
        or any(type(value) is not arb or not value.is_finite() for value in root_box)
    ):
        raise ValueError("regular root box must contain four finite Arb intervals")
    if root_box[3].lower() < 0:
        raise ValueError("regular root multiplier box must be non-negative")

    center = tuple(arb(value.mid()) for value in root_box)
    point_box: BallVector3 = (root_box[0], root_box[1], root_box[2])
    point_center: BallVector3 = (center[0], center[1], center[2])
    multiplier_center = center[_VECTOR_DIMENSION]
    box_evaluation = _evaluate_regular_contact(
        frame, cone, point_box, root_box[_VECTOR_DIMENSION]
    )
    center_evaluation = _evaluate_regular_contact(
        frame, cone, point_center, multiplier_center
    )
    values = (center_evaluation.ellipsoid,) + center_evaluation.stationarity
    jacobian_box = _regular_contact_jacobian(
        frame, cone, point_box, root_box[_VECTOR_DIMENSION]
    )
    midpoint_jacobian = _regular_contact_jacobian(
        frame, cone, point_center, multiplier_center
    ).mid()
    try:
        inverse_enclosure = midpoint_jacobian.inv()
    except ZeroDivisionError as exc:
        raise IntervalCertificationError(
            "regular root midpoint Jacobian is not invertible"
        ) from exc
    try:
        preconditioner = arb_mat(
            [
                [
                    ball_from_float(float(inverse_enclosure[row, column].mid()))
                    for column in range(_VECTOR_DIMENSION + 1)
                ]
                for row in range(_VECTOR_DIMENSION + 1)
            ]
        )
    except (OverflowError, ValueError) as exc:
        raise IntervalCertificationError(
            "regular root inverse has no finite point preconditioner"
        ) from exc
    if not preconditioner.det().is_finite() or preconditioner.det().contains(0):
        raise IntervalCertificationError(
            "regular root point preconditioner is singular"
        )
    image_matrix = krawczyk_image(
        center, values, jacobian_box, root_box, preconditioner
    )
    domain_matrix = arb_mat([[value] for value in root_box])
    if not strictly_contains_vector(domain_matrix, image_matrix):
        raise IntervalCertificationError(
            "regular root Krawczyk image is not strictly inside its box"
        )
    image = tuple(image_matrix[index, 0] for index in range(_VECTOR_DIMENSION + 1))
    return _RegularRootCertificate(
        frame,
        cone,
        root_box,
        image,  # type: ignore[arg-type]
        box_evaluation.cone,
    )


def _reduce_regular_residuals(
    cone: _ConeSection,
    candidates: tuple[_RegularRootCertificate, ...],
) -> _ResidualReduction:
    """Order residuals within one supplied regular-candidate set.

    This helper deliberately has no global or completeness parameter. It only
    reduces candidates whose validity has already been certified. A future
    active-set enumerator must prove that its union covers the whole constrained
    problem before interpreting this local reduction as the global reach.

    Args:
        cone: The cone whose positive cosine normalizes the residual.
        candidates: Non-empty immutable tuple of certified regular candidates.

    Returns:
        Ordered minimum enclosure, normalized negative residual, certified sign,
        and candidate count.

    Raises:
        ValueError: If the candidate container is malformed.
        IntervalCertificationError: If residual intervals overlap in a way that
            prevents one ordered minimum or the normalized sign is unresolved.
    """
    cone.validate()
    if (
        type(candidates) is not tuple
        or not candidates
        or any(
            type(candidate) is not _RegularRootCertificate for candidate in candidates
        )
    ):
        raise ValueError("regular candidates must be a non-empty certificate tuple")
    frame = candidates[0].frame
    if any(
        candidate.frame != frame or candidate.cone != cone for candidate in candidates
    ):
        raise ValueError("regular candidates must belong to one frame and cone")
    residuals = [candidate.cone_residual for candidate in candidates]
    if any(not residual.is_finite() for residual in residuals):
        raise ValueError("regular candidate residuals must be finite intervals")
    minimum_index: int | None = None
    for index, residual in enumerate(residuals):
        if all(
            index == other or residual.upper() < competing.lower()
            for other, competing in enumerate(residuals)
        ):
            minimum_index = index
            break
    if len(residuals) == 1:
        minimum_index = 0
    if minimum_index is None:
        raise IntervalCertificationError(
            "regular candidates do not determine one minimum enclosure"
        )
    minimum = residuals[minimum_index]
    reach = -minimum / ball_from_float(cone.cosine)
    if reach > 0:
        sign = 1
    elif reach < 0:
        sign = -1
    elif reach.is_exact() and reach.is_zero():
        sign = 0
    else:
        raise IntervalCertificationError("regular residual sign is not certified")
    return _ResidualReduction(minimum, reach, sign, len(candidates))


def _fraction_ball(value: Fraction) -> Ball:
    """Convert one exact rational payload field to an Arb point."""
    return arb(f"{value.numerator}/{value.denominator}")


def _dual_vector(witness: DualWitness) -> BallVector3:
    """Construct a dual vector from one exact orthogonality identity."""
    if type(witness) not in {_ProjectedDualWitness, _CrossDualWitness}:
        raise ValueError("dual witness has the wrong private type")
    if (
        type(witness.raw_span) is not tuple
        or len(witness.raw_span) != _VECTOR_DIMENSION
        or type(witness.seed) is not tuple
        or len(witness.seed) != _VECTOR_DIMENSION
    ):
        raise ValueError("dual witness vectors must be exact three-component tuples")
    values = (*witness.raw_span, *witness.seed, witness.scale, witness.multiplier)
    if any(type(value) is not Fraction for value in values):
        raise ValueError("dual witness fields must be exact Fractions")
    if witness.multiplier < 0:
        raise ValueError("dual multiplier must be non-negative")
    span_norm = sum((value * value for value in witness.raw_span), Fraction(0))
    if span_norm <= 0:
        raise ValueError("dual raw span must have positive exact norm")
    if type(witness) is _ProjectedDualWitness:
        seed_dot = sum(
            (witness.seed[index] * witness.raw_span[index] for index in range(3)),
            Fraction(0),
        )
        exact_vector = tuple(
            witness.scale
            * (span_norm * witness.seed[index] - seed_dot * witness.raw_span[index])
            for index in range(3)
        )
    else:
        cross = (
            witness.raw_span[1] * witness.seed[2]
            - witness.raw_span[2] * witness.seed[1],
            witness.raw_span[2] * witness.seed[0]
            - witness.raw_span[0] * witness.seed[2],
            witness.raw_span[0] * witness.seed[1]
            - witness.raw_span[1] * witness.seed[0],
        )
        exact_vector = tuple(witness.scale * value for value in cross)
    if (
        sum(
            (exact_vector[index] * witness.raw_span[index] for index in range(3)),
            Fraction(0),
        )
        != 0
    ):
        raise IntervalCertificationError("dual orthogonality identity did not cancel")
    return tuple(_fraction_ball(value) for value in exact_vector)  # type: ignore[return-value]


def _evaluate_global_reach(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    primal_point: BallVector3 | None,
    dual_witness: DualWitness | None,
    equality_proof: _ContactEqualityProof | None = None,
) -> _GlobalReachProof:
    """Evaluate strict global bounds and an optional exact contact identity.

    A strict primal upper bound proves reach without a dual witness; a strict
    dual lower bound proves miss without a primal witness. ``CONTACT`` requires
    the separate typed equality payload and never follows from a zero-containing
    ordinary bracket. Public eclipse consumers remain disconnected.
    """
    lower: Ball | None = None
    upper: Ball | None = None
    reasons: list[str] = []
    try:
        frame.validate()
        cone.validate()
        if primal_point is None and dual_witness is None and equality_proof is None:
            raise ValueError("global reach needs a primal, dual, or equality witness")
        if equality_proof is not None:
            if type(equality_proof) is not _ContactEqualityProof:
                raise ValueError("contact equality proof has the wrong private type")
            if primal_point is not None:
                raise ValueError("zero-angle line equality has no primal point payload")
            if dual_witness is not None and dual_witness != equality_proof.dual_witness:
                raise ValueError("contact equality dual witness is inconsistent")
            dual_witness = equality_proof.dual_witness
        axis = frame.certified_axis_line()
        metric = frame.metric_ball_matrix()
        cosine = ball_from_float(cone.cosine)
        sine_squared = 1 - cosine * cosine
        if not sine_squared.is_finite() or not (
            sine_squared > 0 or sine_squared.is_exact() and sine_squared.is_zero()
        ):
            raise IntervalCertificationError("cone sine is unresolved")
        sine = sine_squared.sqrt()
        tangent = sine / cosine
        branch = arb(cone.branch_sign)
        radius = ball_from_float(cone.radius_km)

        if primal_point is not None:
            if (
                type(primal_point) is not tuple
                or len(primal_point) != 3
                or any(
                    type(value) is not arb
                    or not value.is_exact()
                    or not value.is_finite()
                    for value in primal_point
                )
            ):
                raise ValueError(
                    "primal point must contain three finite exact Arb points"
                )
            metric_point = _metric_vector(metric, primal_point)
            ellipsoid = _dot(primal_point, metric_point) - 1
            axial = _dot(primal_point, axis.direction)
            relative = tuple(
                primal_point[index] - axis.anchor_km[index] for index in range(3)
            )
            projected = _dot(relative, axis.direction)
            perpendicular = tuple(
                relative[index] - projected * axis.direction[index]
                for index in range(3)
            )
            radial_squared = _dot(perpendicular, perpendicular)
            if not radial_squared.is_finite() or not (
                radial_squared > 0
                or radial_squared.is_exact()
                and radial_squared.is_zero()
            ):
                reasons.append("primal radial distance is unresolved")
            else:
                radial = radial_squared.sqrt()
                nappe_radius = radius + branch * tangent * axial
                primal_objective = cosine * (radial - nappe_radius)
                if ellipsoid <= 0 and nappe_radius >= 0:
                    upper = arb(primal_objective.upper().mid())
                else:
                    reasons.append("primal feasibility is not certified")

        if dual_witness is not None:
            dual_vector = _dual_vector(dual_witness)
            frame_span = tuple(Fraction.from_float(value) for value in frame.axis_span)
            if dual_witness.raw_span != frame_span:
                raise ValueError("dual raw span does not match the contact frame")
            dual_norm_squared = _dot(dual_vector, dual_vector)
            if not dual_norm_squared <= cosine * cosine:
                reasons.append("dual norm is not certified")
            else:
                multiplier = _fraction_ball(dual_witness.multiplier)
                direction = tuple(
                    dual_vector[index]
                    - branch * (sine + multiplier * tangent) * axis.direction[index]
                    for index in range(3)
                )
                try:
                    inverse = metric.inv()
                except ZeroDivisionError as exc:
                    raise IntervalCertificationError(
                        "dual metric inverse failed"
                    ) from exc
                direction_column = arb_mat([[value] for value in direction])
                radicand = (direction_column.transpose() * inverse * direction_column)[
                    0, 0
                ]
                if not radicand.is_finite() or not (
                    radicand > 0 or radicand.is_exact() and radicand.is_zero()
                ):
                    reasons.append("dual support radicand is unresolved")
                else:
                    support = radicand.sqrt()
                    dual_value = (
                        -support
                        - _dot(dual_vector, axis.anchor_km)
                        - cosine * radius
                        - multiplier * radius
                    )
                    lower = arb(dual_value.lower().mid())

        contact = False
        if equality_proof is not None:
            if equality_proof.variant is not _ContactEqualityVariant.ZERO_ANGLE_LINE:
                raise ValueError("contact equality variant is not implemented")
            if not sine_squared.is_exact() or not sine_squared.is_zero():
                raise ValueError("zero-angle contact requires an exact zero cone angle")
            if not radius.is_exact() or not radius.is_zero():
                raise ValueError("zero-angle line contact requires zero cone radius")
            if dual_witness is None:
                raise ValueError("contact equality proof is incomplete")
            dual_vector = _dual_vector(dual_witness)
            if any(not component.is_zero() for component in dual_vector):
                raise ValueError(
                    "zero-angle contact requires an exact zero dual vector"
                )
            if dual_witness.multiplier != 0:
                raise ValueError("zero-angle contact requires a zero dual multiplier")
            metric_direction = _metric_vector(metric, axis.direction)
            metric_anchor = _metric_vector(metric, axis.anchor_km)
            line_a = _dot(axis.direction, metric_direction)
            line_b = _dot(axis.direction, metric_anchor)
            line_c = _dot(axis.anchor_km, metric_anchor) - 1
            discriminant = line_b * line_b - line_a * line_c
            if equality_proof.line_relation is _ZeroAngleLineRelation.CROSSING:
                if discriminant < 0:
                    raise ValueError("zero-angle axis does not reach the ellipsoid")
                if not discriminant > 0:
                    raise IntervalCertificationError(
                        "zero-angle crossing discriminant is unresolved"
                    )
            elif equality_proof.line_relation is _ZeroAngleLineRelation.TANGENT:
                exact_metric = tuple(
                    tuple(Fraction.from_float(value) for value in row)
                    for row in frame.metric_km_minus_2
                )
                exact_point = tuple(
                    Fraction.from_float(value) for value in frame.axis_point_km
                )
                exact_span = tuple(
                    Fraction.from_float(value) for value in frame.axis_span
                )
                exact_a = sum(
                    exact_span[row]
                    * sum(
                        exact_metric[row][column] * exact_span[column]
                        for column in range(_VECTOR_DIMENSION)
                    )
                    for row in range(_VECTOR_DIMENSION)
                )
                exact_b = sum(
                    exact_span[row]
                    * sum(
                        exact_metric[row][column] * exact_point[column]
                        for column in range(_VECTOR_DIMENSION)
                    )
                    for row in range(_VECTOR_DIMENSION)
                )
                exact_c = (
                    sum(
                        exact_point[row]
                        * sum(
                            exact_metric[row][column] * exact_point[column]
                            for column in range(_VECTOR_DIMENSION)
                        )
                        for row in range(_VECTOR_DIMENSION)
                    )
                    - 1
                )
                if exact_b * exact_b - exact_a * exact_c != 0:
                    raise ValueError("zero-angle tangency identity does not hold")
            else:
                raise ValueError("zero-angle line relation has the wrong private type")
            if lower is None or not lower.is_exact() or not lower.is_zero():
                raise ValueError("zero-angle dual equality is not exact")
            contact = True
            lower = arb(0)
            upper = arb(0)

        if contact:
            status = _GlobalReachStatus.CONTACT
            reason = "exact zero-angle line equality proof"
        elif upper is not None and upper < 0:
            status = _GlobalReachStatus.REACH
            reason = "strict primal upper bound"
        elif lower is not None and lower > 0:
            status = _GlobalReachStatus.MISS
            reason = "strict dual lower bound"
        else:
            status = _GlobalReachStatus.UNRESOLVED
            reason = "; ".join(reasons) or "global bounds do not separate a result"
        if lower is not None and upper is not None and lower > upper:
            return _GlobalReachProof(
                _GlobalReachStatus.INVALID,
                lower,
                upper,
                primal_point,
                dual_witness,
                equality_proof,
                "primal and dual bounds are contradictory",
            )
        return _GlobalReachProof(
            status,
            lower,
            upper,
            primal_point,
            dual_witness,
            equality_proof,
            reason,
        )
    except IntervalCertificationError as exc:
        return _GlobalReachProof(
            _GlobalReachStatus.UNRESOLVED,
            lower,
            upper,
            primal_point,
            dual_witness,
            equality_proof,
            str(exc),
        )
    except (ValueError, ZeroDivisionError) as exc:
        return _GlobalReachProof(
            _GlobalReachStatus.INVALID,
            lower,
            upper,
            primal_point,
            dual_witness,
            equality_proof,
            str(exc),
        )


def _validate_vector(vector: object, name: str) -> None:
    """Require one tuple of three finite native Python floats."""
    if type(vector) is not tuple or len(vector) != _VECTOR_DIMENSION:
        raise ValueError(f"{name} must be a tuple of three native floats")
    if any(type(value) is not float or not math.isfinite(value) for value in vector):
        raise ValueError(f"{name} must contain finite native floats")


def _validate_metric(metric: object) -> None:
    """Require an exactly symmetric native-float metric and certify SPD."""
    if type(metric) is not tuple or len(metric) != _VECTOR_DIMENSION:
        raise ValueError("ellipsoid metric must be a tuple of three rows")
    for row in metric:
        _validate_vector(row, "ellipsoid metric row")
    for row in range(_VECTOR_DIMENSION):
        for column in range(row):
            if metric[row][column] != metric[column][row]:
                raise ValueError("ellipsoid metric must be exactly symmetric")
    exact = tuple(tuple(Fraction.from_float(value) for value in row) for row in metric)
    leading_minors = (
        exact[0][0],
        exact[0][0] * exact[1][1] - exact[0][1] * exact[1][0],
        exact[0][0] * (exact[1][1] * exact[2][2] - exact[1][2] * exact[2][1])
        - exact[0][1] * (exact[1][0] * exact[2][2] - exact[1][2] * exact[2][0])
        + exact[0][2] * (exact[1][0] * exact[2][1] - exact[1][1] * exact[2][0]),
    )
    if any(minor <= 0 for minor in leading_minors):
        raise ValueError("ellipsoid metric must be positive definite")
    interval_cholesky(
        arb_mat([[ball_from_float(value) for value in row] for row in metric])
    )


def _ball_vector(vector: Vector3) -> BallVector3:
    """Convert one validated native-float vector to exact Arb inputs."""
    return tuple(ball_from_float(value) for value in vector)  # type: ignore[return-value]
