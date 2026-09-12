# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Private certified tuples for the four secular perturbing planets.

This module is intentionally not connected to a public calculation path.  It
implements the internal B-11 algebra contract: exact binary64 inputs are
promoted to fresh Arb balls, source-frame elements are converted through the
relative ERFA ecliptic frame, and two independent precision passes must select
the same five native floats.
"""

from __future__ import annotations

import hashlib
import math
import struct
import threading
from dataclasses import dataclass
from typing import TypeAlias

import erfa
from flint import arb, ctx

from .constants import JUPITER, NEPTUNE, SATURN, URANUS
from .exceptions import CalculationError
from .intervals import ball_from_float, certified_float
from .planetary_mean_elements import (
    _MeanElementCoefficients,
    _mean_element_coefficients,
)

_J2000 = 2451545.0
_CENTURY_DAYS = 36525.0
_RUNTIME_SCHEDULES = ((192, 384, 768, 1536, 3072, 4096), (256, 512, 1024, 2048, 4096))
_PRECISION_LOCK = threading.RLock()

Vector: TypeAlias = tuple[arb, arb, arb]
Matrix: TypeAlias = tuple[
    tuple[arb, arb, arb], tuple[arb, arb, arb], tuple[arb, arb, arb]
]


@dataclass(frozen=True, slots=True)
class _PerturberProof:
    """Optional internal proof intermediates retained for focused tests."""

    proof_class: str
    evec_date: Vector
    pole_date: Vector
    evec_j2000: Vector
    pole_j2000: Vector
    frame_date: Matrix | None
    frame_j2000: Matrix | None
    frame: Matrix
    inverse_frame: Matrix


def _canonical_coefficients_bytes() -> bytes:
    """Encode the registered rows using the B-11 canonical binary format."""
    chunks: list[bytes] = []
    for body_id in (JUPITER, SATURN, URANUS, NEPTUNE):
        row = _mean_element_coefficients(body_id)
        assert row is not None
        chunks.append(struct.pack(">q", body_id))
        for name in ("a", "e", "i", "node", "perihelion"):
            values = getattr(row, name)
            encoded_name = name.encode("utf-8")
            chunks.append(struct.pack(">I", len(encoded_name)))
            chunks.append(encoded_name)
            chunks.append(struct.pack(">I", len(values)))
            chunks.extend(struct.pack(">d", value) for value in values)
    return b"".join(chunks)


def _coefficient_source_sha256() -> str:
    """Return the digest of the canonical coefficient serialization."""
    return hashlib.sha256(_canonical_coefficients_bytes()).hexdigest()


def _poly_ball(coefficients: tuple[float, ...], t: arb) -> arb:
    """Evaluate an ascending-power row with descending Horner steps."""
    result = arb(0)
    for coefficient in reversed(coefficients):
        result = result * t + ball_from_float(coefficient)
    return result


def _cross(left: Vector, right: Vector) -> Vector:
    """Cross two vectors in the registered left-to-right operation order."""
    return (
        (left[1] * right[2]) - (left[2] * right[1]),
        (left[2] * right[0]) - (left[0] * right[2]),
        (left[0] * right[1]) - (left[1] * right[0]),
    )


def _dot(left: Vector, right: Vector) -> arb:
    """Three-term dot product with the specified grouping."""
    return ((left[0] * right[0]) + (left[1] * right[1])) + (left[2] * right[2])


def _norm(value: Vector) -> arb:
    """Euclidean norm with the specified grouping."""
    return (
        ((value[0] * value[0]) + (value[1] * value[1])) + value[2] * value[2]
    ).sqrt()


def _require_vector(value: Vector, name: str) -> arb:
    """Validate a finite vector norm and return it."""
    norm = _norm(value)
    if not norm.is_finite() or norm <= 0:
        raise CalculationError(f"{name} has no positive finite norm")
    return norm


def _normalize(value: Vector, name: str) -> Vector:
    """Normalize an Arb vector only after proving a positive finite norm."""
    denominator = _require_vector(value, name)
    return (value[0] / denominator, value[1] / denominator, value[2] / denominator)


def _frame_rows(jd_tt: float) -> Matrix:
    """Construct a mean ecliptic frame directly from the two ERFA pole calls."""
    epj = 2000.0 + (jd_tt - 2451545.0) / 365.25
    ecliptic_values = erfa.ltpecl(epj)
    equator_values = erfa.ltpequ(epj)
    ecliptic_pole: Vector = (
        ball_from_float(float(ecliptic_values[0])),
        ball_from_float(float(ecliptic_values[1])),
        ball_from_float(float(ecliptic_values[2])),
    )
    equator_pole: Vector = (
        ball_from_float(float(equator_values[0])),
        ball_from_float(float(equator_values[1])),
        ball_from_float(float(equator_values[2])),
    )
    z = _normalize(ecliptic_pole, "ecliptic pole")
    x = _normalize(_cross(z, equator_pole), "equinox cross product")
    y = _normalize(_cross(z, x), "frame y row")
    rows: Matrix = (x, y, z)
    for first in rows:
        norm = _require_vector(first, "frame row")
        if not (norm * norm).contains(1):
            raise CalculationError("ecliptic frame row is not unit length")
    for row_index in range(3):
        for column_index in range(row_index):
            if not _dot(rows[row_index], rows[column_index]).contains(0):
                raise CalculationError("ecliptic frame rows are not orthogonal")
    if _dot(_cross(x, y), z) <= 0:
        raise CalculationError("ecliptic frame is not right-handed")
    return rows


def _identity() -> Matrix:
    """Return the literal binary64 identity as exact Arb balls."""
    return (
        (arb(1), arb(0), arb(0)),
        (arb(0), arb(1), arb(0)),
        (arb(0), arb(0), arb(1)),
    )


def _transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[row][column] for row in range(3)) for column in range(3))  # type: ignore[return-value]


def _matrix_product(left: Matrix, right: Matrix) -> Matrix:
    """Multiply matrices using the specified grouped three-term dot product."""
    return tuple(
        tuple(
            ((left[row][0] * right[0][column]) + (left[row][1] * right[1][column]))
            + (left[row][2] * right[2][column])
            for column in range(3)
        )
        for row in range(3)
    )  # type: ignore[return-value]


def _matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(_dot(matrix[row], vector) for row in range(3))  # type: ignore[return-value]


def _relative_frames(
    jd_tt: float, t: arb
) -> tuple[Matrix, Matrix, Matrix, Matrix | None, Matrix | None]:
    """Return relative frame and reverse, with literal identity at J2000."""
    if jd_tt == _J2000 and t.is_exact() and t.is_zero():
        identity = _identity()
        return identity, identity, identity, None, None
    frame0 = _frame_rows(_J2000)
    frame_date = _frame_rows(jd_tt)
    relative = _matrix_product(frame0, _transpose(frame_date))
    reverse = _matrix_product(frame_date, _transpose(frame0))
    return relative, reverse, relative, frame_date, frame0


def _source_vectors(
    row: _MeanElementCoefficients, t: arb
) -> tuple[arb, Vector, Vector]:
    """Evaluate full elements and make physical eccentricity and pole vectors."""
    a = _poly_ball(row.a, t)
    eccentricity = _poly_ball(row.e, t)
    inclination = _poly_ball(row.i, t) * (arb.pi() / 180)
    node = _poly_ball(row.node, t) * (arb.pi() / 180)
    perihelion = (_poly_ball(row.perihelion, t) - _poly_ball(row.node, t)) * (
        arb.pi() / 180
    )
    cos_node, sin_node = node.cos(), node.sin()
    cos_i, sin_i = inclination.cos(), inclination.sin()
    cos_perihelion, sin_perihelion = perihelion.cos(), perihelion.sin()
    p: Vector = (
        cos_node * cos_perihelion - sin_node * cos_i * sin_perihelion,
        sin_node * cos_perihelion + cos_node * cos_i * sin_perihelion,
        sin_i * sin_perihelion,
    )
    pole = (sin_i * sin_node, -sin_i * cos_node, cos_i)
    return a, (eccentricity * p[0], eccentricity * p[1], eccentricity * p[2]), pole


def _absolute_ball(value: arb) -> arb:
    if value > 0:
        return value
    if value < 0:
        return -value
    return value.union(-value)


def _proved_sign(value: arb, name: str) -> int:
    """Return a proved sign, distinguishing exact zero from an unresolved ball."""
    if value.is_exact() and value.is_zero():
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise CalculationError(f"{name} sign is not separated from zero")


def _angle_chart(pair_x: arb, pair_y: arb) -> arb:
    """Evaluate a finite local atan2 chart or fail closed."""
    x_sign = _proved_sign(pair_x, "atan2 x")
    y_sign = _proved_sign(pair_y, "atan2 y")
    if x_sign == 0 and y_sign == 0:
        raise CalculationError("angle pair is an unresolved origin")
    if x_sign == 0:
        return arb.pi() / 2 if y_sign > 0 else -arb.pi() / 2
    if y_sign == 0:
        return arb(0) if x_sign > 0 else arb.pi()

    abs_x = _absolute_ball(pair_x)
    abs_y = _absolute_ball(pair_y)
    difference = abs_x - abs_y
    if difference > 0:
        denominator = pair_x
        ratio = _absolute_ball(pair_y / denominator)
        minor = ratio.atan()
        if x_sign > 0:
            return minor if y_sign > 0 else -minor
        return arb.pi() - minor if y_sign > 0 else -arb.pi() + minor
    if difference < 0:
        denominator = pair_y
        ratio = _absolute_ball(pair_x / denominator)
        minor = ratio.atan()
        if y_sign > 0:
            return arb.pi() / 2 - minor if x_sign > 0 else arb.pi() / 2 + minor
        return -arb.pi() / 2 + minor if x_sign > 0 else -arb.pi() / 2 - minor

    # At a seam both denominators are nonzero.  Intersect, never hull, the two
    # valid charts evaluated against this same correlated pair.
    # Evaluate each
    # formula explicitly at the seam where |x| == |y|.
    horizontal_minor = (_absolute_ball(pair_y / pair_x)).atan()
    vertical_minor = (_absolute_ball(pair_x / pair_y)).atan()
    first = (
        horizontal_minor
        if x_sign > 0 and y_sign > 0
        else -horizontal_minor
        if x_sign > 0
        else arb.pi() - horizontal_minor
        if y_sign > 0
        else -arb.pi() + horizontal_minor
    )
    second = (
        arb.pi() / 2 - vertical_minor
        if y_sign > 0 and x_sign > 0
        else arb.pi() / 2 + vertical_minor
        if y_sign > 0
        else -arb.pi() / 2 + vertical_minor
        if x_sign > 0
        else -arb.pi() / 2 - vertical_minor
    )
    try:
        result = first.intersection(second)
    except ValueError as exc:
        raise CalculationError("angle chart intersection is empty") from exc
    if not result.is_finite() or result.lower() > result.upper():
        raise CalculationError("angle chart intersection is invalid")
    return result


def _axis_angle_degrees(pair_x: arb, pair_y: arb) -> float | None:
    """Return a canonical degree only for a proved exact axis pair."""
    x_sign = _proved_sign(pair_x, "atan2 x")
    y_sign = _proved_sign(pair_y, "atan2 y")
    if x_sign == 0 and y_sign == 0:
        raise CalculationError("angle pair is an undefined exact origin")
    if x_sign == 0:
        return 90.0 if y_sign > 0 else 270.0
    if y_sign == 0:
        return 0.0 if x_sign > 0 else 180.0
    return None


def _angle_degrees(pair_x: arb, pair_y: arb) -> float:
    """Certify an ``atan2(pair_y, pair_x)`` angle in canonical degrees."""
    axis = _axis_angle_degrees(pair_x, pair_y)
    if axis is not None:
        return axis
    return _normalized_degrees(_angle_chart(pair_x, pair_y))


def _normalized_degrees(angle: arb) -> float:
    """Certify an angle in degrees after one exact integer-turn translation."""
    midpoint = float(angle.mid()) * (180.0 / math.pi)
    if not math.isfinite(midpoint):
        raise CalculationError("angle is not finite")
    turns = math.floor(midpoint / 360.0)
    normalized = angle * (180.0 / math.pi) - 360.0 * turns
    candidate = float(normalized.mid())
    candidate %= 360.0
    # Keep the interval near the candidate's normalized branch.  A value at the
    # upper endpoint is represented by zero, not by 360.
    if candidate == 0.0 and normalized.mid() > 359:
        normalized = normalized - 360.0
    try:
        return float(certified_float(normalized))
    except Exception as exc:
        raise CalculationError("angle does not select one binary64 cell") from exc


def _recover_tuple(
    a: arb, evec: Vector, pole_value: Vector
) -> tuple[tuple[float, float, float, float, float], str, Vector, Vector]:
    """Recover the target tuple and canonical exact-degeneracy proof class."""
    pole_norm = _require_vector(pole_value, "transformed pole")
    pole: Vector = (
        pole_value[0] / pole_norm,
        pole_value[1] / pole_norm,
        pole_value[2] / pole_norm,
    )
    eccentricity = _norm(evec)
    if not eccentricity.is_finite() or eccentricity < 0:
        raise CalculationError("eccentricity vector is not finite")

    pole_z = pole[2]
    if pole_z < -1 or pole_z > 1:
        raise CalculationError("pole is outside the acos domain")
    inclination = pole_z.acos()
    sin_i = (pole[0] * pole[0] + pole[1] * pole[1]).sqrt()
    exact_zero_e = eccentricity.is_exact() and eccentricity.is_zero()
    exact_zero_i = sin_i.is_exact() and sin_i.is_zero()
    exact_retrograde = exact_zero_i and pole[2].is_exact() and pole[2] == -1
    if exact_zero_i:
        if exact_zero_e:
            omega_deg = 0.0
        elif exact_retrograde:
            omega_deg = _angle_degrees(evec[0], -evec[1])
        else:
            omega_deg = _angle_degrees(evec[0], evec[1])
        tuple_value = (
            float(certified_float(a)),
            float(certified_float(eccentricity)),
            180.0 if exact_retrograde else 0.0,
            0.0,
            omega_deg,
        )
        return (
            tuple_value,
            "canonical-pole" if exact_retrograde else "canonical-equator",
            pole,
            evec,
        )

    omega_node_deg = _angle_degrees(-pole[1], pole[0])
    n_sin = sin_i
    node_basis: Vector = (-pole[1] / n_sin, pole[0] / n_sin, arb(0))
    quadrature = _cross(pole, node_basis)
    if exact_zero_e:
        omega_deg = 0.0
    else:
        perihelion: Vector = (
            evec[0] / eccentricity,
            evec[1] / eccentricity,
            evec[2] / eccentricity,
        )
        omega_deg = _angle_degrees(
            _dot(perihelion, node_basis), _dot(perihelion, quadrature)
        )
    result = (
        float(certified_float(a)),
        float(certified_float(eccentricity)),
        float(certified_float(inclination * (180.0 / math.pi))),
        omega_node_deg,
        omega_deg,
    )
    return result, "zero-eccentricity" if exact_zero_e else "ordinary", pole, evec


def _one_precision_pass(
    body_id: int, jd_tt: float, bits: int
) -> tuple[tuple[float, float, float, float, float], _PerturberProof]:
    """Run one fully fresh Arb/ERFA candidate selection at ``bits`` precision."""
    with _precision(bits):
        row = _mean_element_coefficients(body_id)
        if row is None:
            raise CalculationError("body is not a registered perturbing planet")
        if not math.isfinite(jd_tt):
            raise CalculationError("Julian Day must be finite")
        t = (ball_from_float(jd_tt) - ball_from_float(_J2000)) / ball_from_float(
            _CENTURY_DAYS
        )
        a, evec_date, pole_date = _source_vectors(row, t)
        relative, reverse, frame, frame_date, frame_j2000 = _relative_frames(jd_tt, t)
        evec_j2000 = _matrix_vector(relative, evec_date)
        pole_j2000 = _matrix_vector(relative, pole_date)
        result, proof_class, pole, evec = _recover_tuple(a, evec_j2000, pole_j2000)
        proof = _PerturberProof(
            proof_class=proof_class,
            evec_date=evec_date,
            pole_date=pole_date,
            evec_j2000=evec,
            pole_j2000=pole,
            frame_date=frame_date,
            frame_j2000=frame_j2000,
            frame=frame,
            inverse_frame=reverse,
        )
        return result, proof


class _precision:
    """Save and restore the process-global Arb precision even on failure."""

    def __init__(self, bits: int) -> None:
        self.bits = bits
        self.previous: int | None = None

    def __enter__(self) -> None:
        self.previous = ctx.prec
        ctx.prec = self.bits

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        assert self.previous is not None
        ctx.prec = self.previous


def _certified_perturber_tuple(
    body_id: int, jd_tt: float, *, return_proof: bool = False
) -> (
    tuple[float, float, float, float, float]
    | tuple[tuple[float, float, float, float, float], _PerturberProof]
):
    """Return the certified anchored-J2000 tuple for one perturber.

    This private evaluator is initially unconnected to all public consumers.
    Both complete precision schedules are independent and serialized because
    ``flint.ctx.prec`` is process-global.
    """
    if not isinstance(body_id, int) or isinstance(body_id, bool):
        raise CalculationError("body identifier must be an integer")
    if not math.isfinite(jd_tt):
        raise CalculationError("Julian Day must be finite")
    with _PRECISION_LOCK:
        successes: list[
            tuple[tuple[float, float, float, float, float], _PerturberProof]
        ] = []
        for schedule in _RUNTIME_SCHEDULES:
            selected: (
                tuple[tuple[float, float, float, float, float], _PerturberProof] | None
            ) = None
            for bits in schedule:
                try:
                    selected = _one_precision_pass(body_id, jd_tt, bits)
                    break
                except (CalculationError, ValueError, ArithmeticError):
                    continue
            if selected is None:
                raise CalculationError("perturber tuple was not certified by 4096 bits")
            successes.append(selected)
        first, second = successes
        if first[0] != second[0] or first[1].proof_class != second[1].proof_class:
            raise CalculationError("independent perturber passes disagree")
        return first if return_proof else first[0]


__all__ = ["_certified_perturber_tuple"]
