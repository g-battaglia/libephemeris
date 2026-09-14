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
from typing import Final

from flint import arb, arb_mat

from .intervals import (
    Ball,
    BallMatrix,
    IntervalCertificationError,
    ball_from_float,
    interval_cholesky,
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


@dataclass(frozen=True, slots=True)
class _RegularContactEvaluation:
    """Regular cone/ellipsoid residuals and stationarity vector."""

    ellipsoid: Ball
    cone: Ball
    nappe_radius: Ball
    radial_distance: Ball
    stationarity: BallVector3


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


def _dot(left: tuple[Ball, ...], right: tuple[Ball, ...]) -> Ball:
    """Return one outward-rounded Euclidean dot product."""
    return sum((left[index] * right[index] for index in range(len(left))), arb(0))


def _metric_vector(metric: BallMatrix, vector: BallVector3) -> BallVector3:
    """Multiply one three-dimensional Arb matrix/vector pair."""
    column = arb_mat([[value] for value in vector])
    product = metric * column
    return tuple(product[index, 0] for index in range(_VECTOR_DIMENSION))  # type: ignore[return-value]


def _evaluate_regular_contact(
    frame: _EllipsoidContactFrame,
    cone: _ConeSection,
    point_km: BallVector3,
    ellipsoid_multiplier: Ball,
) -> _RegularContactEvaluation:
    """Evaluate equations (1), (2), and regular-contact stationarity.

    This evaluator covers only ``rho > 0`` and ``r > 0``. It does not solve the
    equations, classify a global minimum, or handle the radial subgradient and
    active nappe boundary. A caller must prove the two strict domain conditions
    before applying Krawczyk to this smooth system.

    Args:
        frame: Validated private ellipsoid and axis source frame.
        cone: One validated physical cone section.
        point_km: Three-dimensional interval point.
        ellipsoid_multiplier: Non-negative KKT multiplier enclosure.

    Returns:
        Outward enclosures of the ellipsoid equation, cone equation, nappe
        radius, radial distance, and three stationarity components.

    Raises:
        ValueError: If a point/multiplier/cone input is malformed or non-finite.
        IntervalCertificationError: If ``rho`` or the nappe radius is not proved
            strictly positive.
    """
    frame.validate()
    cone.validate()
    if (
        type(point_km) is not tuple
        or len(point_km) != _VECTOR_DIMENSION
        or any(not value.is_finite() for value in point_km)
    ):
        raise ValueError("contact point must contain three finite Arb intervals")
    if not ellipsoid_multiplier.is_finite() or ellipsoid_multiplier < 0:
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
    if not nappe_radius > 0:
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
    interval_cholesky(
        arb_mat([[ball_from_float(value) for value in row] for row in metric])
    )


def _ball_vector(vector: Vector3) -> BallVector3:
    """Convert one validated native-float vector to exact Arb inputs."""
    return tuple(ball_from_float(value) for value in vector)  # type: ignore[return-value]
